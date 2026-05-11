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


#ifndef _ICHIMAIN_H_
#define _ICHIMAIN_H_

#include "strutil.h"
#include "ichicomn.h"

#define ESC_KEY       27

#define INCHI_SEGM_BUFLEN  524288
#define PRINT_INCHI_MAX_TAG_LEN 64
typedef struct tagLine
{
    char *str;
    int   len;
    int   len_alloc;
    int   c;
} SEGM_LINE;
/* for DisplayTheWholeStructure() */
#define COMP_ORIG_0_MAIN  0x0001
#define COMP_ORIG_0_RECN  0x0002
#define COMP_PREP_0_MAIN  0x0004
#define COMP_PREP_0_RECN  0x0008
#define COMP_ORIG_1_MAIN  0x0010
#define COMP_ORIG_1_RECN  0x0020


/* STRUCT_DATA */


typedef struct tagStructData
{
    unsigned long ulStructTime;
    int           nErrorCode;
    int           nErrorType;
    int           nStructReadError;
    char          pStrErrStruct[STR_ERR_LEN];
    long          fPtrStart;  /* or number of processed structures */
    long          fPtrEnd;    /* or number of errors */
    int           bUserQuit;
    int           bUserQuitComponent;
    int           bUserQuitComponentDisplay;
    int           bChiralFlag;

    /* information related to normal or disconnected layers */
    int           num_taut[INCHI_NUM];
    int           num_non_taut[INCHI_NUM];
    INCHI_MODE     bTautFlags[INCHI_NUM];        /* reconnected does not have TG_FLAG_DISCONNECT_COORD_DONE flag */
    INCHI_MODE     bTautFlagsDone[INCHI_NUM];    /* reconnected does not have TG_FLAG_DISCONNECT_COORD_DONE flag */
    int           num_components[INCHI_NUM];    /* number of allocated INChI, INChI_Aux data structures */
    /* debugging info */
#if ( bRELEASE_VERSION == 0 )
    int           bExtract;
#endif
} STRUCT_DATA;


/* Convenience storage for InChI serialization control data */
typedef struct tagINCHI_OUT_CTL
{
    int ATOM_MODE;
    int TAUT_MODE;

    int *pSortPrintINChIFlags;

    int bOverflow;
    int bAlways;
    int bOutputType;
    int bOutType;
    int bPlainTextTags;
    int bOmitRepetitions;
    int bUseMulipliers;
    int bNonTautNonIsoIdentifierNotEmpty;
    int bNonTautIsoIdentifierNotEmpty;
    int bSecondNonTautPass;
    int bTautomericOutputAllowed;
    int bTautomeric;
    int bNonTautomeric;
    int bNonTautIsIdenticalToTaut;
    int bFhTag;
    int bRelRac;
    int bAbcNumbers;
    int bIsotopic;
    int bPolymers;

    int iCurTautMode;

    int num_components;
    int  nNumRemovedProtons;
    int nTag;
    int bTag1;
    int bTag2;
    int bTag3;
    int tot_len;
    int tot_len2;

    int nCurINChISegment;
    int nSegmAction;

    int num_comp[TAUT_NUM];
    int num_iso_H[NUM_H_ISOTOPES];
    int bAtomEqu[TAUT_NUM];
    int bTautEqu[TAUT_NUM];
    int bInvStereo[TAUT_NUM];
    int bInvStereoOrigNumb[TAUT_NUM];
    int bRacemicStereo[TAUT_NUM];
    int bRelativeStereo[TAUT_NUM];
    int bIsotopicOrigNumb[TAUT_NUM];
    int bIsotopicAtomEqu[TAUT_NUM];
    int bIsotopicTautEqu[TAUT_NUM];
    int bInvIsotopicStereo[TAUT_NUM];
    int bInvIsotopicStereoOrigNumb[TAUT_NUM];
    int bIsotopicRacemicStereo[TAUT_NUM];
    int bIsotopicRelativeStereo[TAUT_NUM];
    int bIgn_UU_Sp3[TAUT_NUM];
    int bIgn_UU_Sp2[TAUT_NUM];
    int bIgn_UU_Sp3_Iso[TAUT_NUM];
    int bIgn_UU_Sp2_Iso[TAUT_NUM];
    int bChargesRadVal[TAUT_NUM];
    int bOrigCoord[TAUT_NUM];

    char sDifSegs[DIFL_LENGTH][DIFS_LENGTH];
    char szTag1[PRINT_INCHI_MAX_TAG_LEN];
    char szTag2[PRINT_INCHI_MAX_TAG_LEN];
    char szTag3[PRINT_INCHI_MAX_TAG_LEN];

    int n_pzz; 
    int n_zy;

    INCHI_SORT   **pINChISortTautAndNonTaut;
    INCHI_SORT   *pINChISort;
    INCHI_SORT   *pINChISort2;
}
INCHI_OUT_CTL;


/* Context of ProcessOneStructureEx used in calls to ProcessOneStructureExCore */
typedef struct tagPOSEContext
{
    STRUCT_DATA sd;
    INPUT_PARMS ip;
    char szTitle[MAX_SDF_HEADER + MAX_SDF_VALUE + 256];
    PINChI2 *pINChI2[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux2[INCHI_NUM];
    INCHI_IOSTREAM *inp_file;
    INCHI_IOSTREAM inchi_file[3];
    INCHI_IOSTREAM *log_file;
    INCHI_IOSTREAM *out_file;
    INCHI_IOSTREAM *prb_file;
    ORIG_ATOM_DATA OrigAtData;
    ORIG_ATOM_DATA *orig_inp_data;
    ORIG_ATOM_DATA PrepAtData[2];
    ORIG_ATOM_DATA *prep_inp_data;
    long num_inp;
    INCHI_IOS_STRING temp_string_container;
    INCHI_IOS_STRING *strbuf;
    unsigned char save_opt_bits;
} POSEContext;
int  POSEContext_Init(POSEContext *context,
                    STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                    PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                    INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                    INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file,
                    ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                    long num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits);
void POSEContext_Free(POSEContext *context);
void POSEContext_DebugPrint(POSEContext *context);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/* MAIN_LOOP_ACTION */

typedef enum MAIN_LOOP_ACTION
{
    DO_NEXT_STEP,
    DO_BREAK_MAIN_LOOP,
    DO_EXIT_FUNCTION,
    DO_CONTINUE_MAIN_LOOP
}
MAIN_LOOP_ACTION;


/* GENERAL PROCESSING STEPS */


int ProcessSingleInputFile( int argc, char *argv[] );
int ProcessMultipleInputFiles( int argc, char *argv[] );
int ReadCommandLineParms( int argc, const char *argv[], INPUT_PARMS *ip,
                          char *szSdfDataValue, unsigned long *ulDisplTime,
                          int bReleaseVersion, INCHI_IOSTREAM *log_file );
void HelpCommandLineParms( INCHI_IOSTREAM *f );
int OpenFiles( FILE **inp_file, FILE **out_file, FILE **log_file, FILE **prb_file, INPUT_PARMS *ip );
int PrintInputParms( INCHI_IOSTREAM *log_file, INPUT_PARMS *ip );
int SortAndPrintINChI( struct tagCANON_GLOBALS *pCG,
                       INCHI_IOSTREAM *out_file,
                       INCHI_IOS_STRING *strbuf,
                       INCHI_IOSTREAM *log_file,
                       INPUT_PARMS *ip,
                       ORIG_ATOM_DATA *orig_inp_data,
                       ORIG_ATOM_DATA *prep_inp_data,
                       COMP_ATOM_DATA composite_norm_data[INCHI_NUM][TAUT_NUM + 1],
                       ORIG_STRUCT *pOrigStruct,
                       int num_components[INCHI_NUM],
                       int num_non_taut[INCHI_NUM],
                       int num_taut[INCHI_NUM],
                       INCHI_MODE bTautFlags[INCHI_NUM],
                       INCHI_MODE bTautFlagsDone[INCHI_NUM],
                       NORM_CANON_FLAGS *pncFlags,
                       long num_inp,
                       PINChI2 *pINChI[INCHI_NUM],
                       PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                       int *pSortPrintINChIFlags,
                       unsigned char save_opt_bits );
void FreeAllINChIArrays( PINChI2 *pINChI[INCHI_NUM],
                         PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                         int num_components[2] );
void FreeINChIArrays( PINChI2 *pINChI,
                      PINChI_Aux2 *pINChI_Aux,
                      int num_components );
int ReadTheStructure( struct tagINCHI_CLOCK *ic,
                      STRUCT_DATA *sd,
                      INPUT_PARMS *ip,
                      INCHI_IOSTREAM *inp_file,
                      ORIG_ATOM_DATA *orig_inp_data,
                      int inp_index,
                      int *out_index );
int TreatErrorsInReadTheStructure( STRUCT_DATA *sd,
                                   INPUT_PARMS *ip,
                                   int nLogMask,
                                   INCHI_IOSTREAM *inp_file,
                                   INCHI_IOSTREAM *log_file,
                                   INCHI_IOSTREAM *out_file,
                                   INCHI_IOSTREAM *prb_file,
                                   ORIG_ATOM_DATA *orig_inp_data,
                                   long *num_inp );
int GetOneComponent( struct tagINCHI_CLOCK *ic,
                     STRUCT_DATA *sd,
                     INPUT_PARMS *ip,
                     INCHI_IOSTREAM *log_file,
                     INCHI_IOSTREAM *out_file,
                     INP_ATOM_DATA *inp_cur_data,
                     ORIG_ATOM_DATA *orig_inp_data,
                     int i, long num_inp );
int CreateOneComponentINChI( struct tagCANON_GLOBALS *pCG,
                             struct tagINCHI_CLOCK *ic,
                             STRUCT_DATA *sd,
                             INPUT_PARMS *ip,
                             INP_ATOM_DATA *inp_cur_data,
                             ORIG_ATOM_DATA *orig_inp_data,
                             PINChI2 *pINChI,
                             PINChI_Aux2 *pINChI_Aux,
                             int iINChI,
                             int i, long num_inp,
                             INP_ATOM_DATA **inp_norm_data,
                             NORM_CANON_FLAGS *pncFlags,
                             INCHI_IOSTREAM *log_file );
int TreatErrorsInCreateOneComponentINChI( STRUCT_DATA *sd,
                                          INPUT_PARMS *ip,
                                          ORIG_ATOM_DATA *orig_inp_data,
                                          int i,
                                          long num_inp,
                                          INCHI_IOSTREAM *inp_file,
                                          INCHI_IOSTREAM *log_file,
                                          INCHI_IOSTREAM *out_file,
                                          INCHI_IOSTREAM *prb_file );
int TreatCreateINChIWarning( STRUCT_DATA *sd,
                             INPUT_PARMS *ip,
                             ORIG_ATOM_DATA *orig_inp_data,
                             long num_inp,
                             INCHI_IOSTREAM *inp_file,
                             INCHI_IOSTREAM *log_file,
                             INCHI_IOSTREAM *out_file,
                             INCHI_IOSTREAM *prb_file );
int GetProcessingWarningsOneComponentInChI( INChI *cur_INChI[],
                                            INP_ATOM_DATA **inp_norm_data,
                                            STRUCT_DATA *sd,
                                            int bNoWarnings );
int GetTheNextRecordOfInputFile( struct tagINCHI_CLOCK *ic,
                                 STRUCT_DATA *sd, INPUT_PARMS *ip,
                                 char *szTitle,
                                 INCHI_IOSTREAM *inp_file,
                                 INCHI_IOSTREAM *plog,
                                 INCHI_IOSTREAM *pout,
                                 INCHI_IOSTREAM *pprb,
                                 ORIG_ATOM_DATA *orig_inp_data,
                                 long *num_inp,
                                 STRUCT_FPTRS *pStructPtrs,
                                 int *nRet,
                                 int *have_err_in_GetOneStructure,
                                 long *num_err,
                                 int output_error_inchi );
int CalcAndPrintINCHIAndINCHIKEY( struct tagINCHI_CLOCK *ic,
                                  CANON_GLOBALS *CG,
                                  STRUCT_DATA *sd,
                                  INPUT_PARMS *ip,
                                  char *szTitle,
                                  PINChI2 *pINChI[INCHI_NUM],
                                  PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                                  INCHI_IOSTREAM *inp_file,
                                  INCHI_IOSTREAM *plog,
                                  INCHI_IOSTREAM *pout,
                                  INCHI_IOSTREAM *pprb,
                                  ORIG_ATOM_DATA *orig_inp_data,
                                  ORIG_ATOM_DATA *prep_inp_data,
                                  long *num_inp,
                                  STRUCT_FPTRS *pStructPtrs,
                                  int *nRet,
                                  int have_err_in_GetOneStructure,
                                  long *num_err,
                                  int output_error_inchi,
                                  INCHI_IOS_STRING *strbuf,
                                  unsigned long *pulTotalProcessingTime,
                                  char *pLF, char *pTAB,
                                  char *ikey, int silent );
int GetOneStructure( struct tagINCHI_CLOCK *ic,
                     STRUCT_DATA *sd,
                     INPUT_PARMS *ip,
                     char *szTitle,
                     INCHI_IOSTREAM *inp_file,
                     INCHI_IOSTREAM *log_file,
                     INCHI_IOSTREAM *out_file,
                     INCHI_IOSTREAM *prb_file,
                     ORIG_ATOM_DATA *orig_inp_data,
                     long *num_inp,
                     STRUCT_FPTRS *struct_fptrs );
int ProcessOneStructure( struct tagINCHI_CLOCK *ic,
                         struct tagCANON_GLOBALS *pCG,
                         STRUCT_DATA *sd,
                         INPUT_PARMS *ip,
                         char *szTitle,
                         PINChI2 *pINChI2[INCHI_NUM],
                         PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                         INCHI_IOSTREAM *inp_file,
                         INCHI_IOSTREAM *log_file,
                         INCHI_IOSTREAM *out_file,
                         INCHI_IOSTREAM *prb_file,
                         ORIG_ATOM_DATA *orig_inp_data,
                         ORIG_ATOM_DATA *prep_inp_data,
                         long num_inp,
                         INCHI_IOS_STRING *strbuf,
                         unsigned char save_opt_bits );
int ProcessOneStructureEx( struct tagINCHI_CLOCK *ic,
                           struct tagCANON_GLOBALS *pCG,
                           STRUCT_DATA *sd,
                           INPUT_PARMS *ip,
                           char *szTitle,
                           PINChI2 *pINChI2[INCHI_NUM],
                           PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                           INCHI_IOSTREAM *inp_file,
                           INCHI_IOSTREAM *log_file,
                           INCHI_IOSTREAM *out_file,
                           INCHI_IOSTREAM *prb_file,
                           ORIG_ATOM_DATA *orig_inp_data,
                           ORIG_ATOM_DATA *prep_inp_data,
                           long num_inp,
                           INCHI_IOS_STRING *strbuf,
                           unsigned char save_opt_bits );
int PreprocessPolymerCRUData( struct tagINCHI_CLOCK    *ic,
                              struct tagCANON_GLOBALS  *CG,
                              STRUCT_DATA              *sd,
                              INPUT_PARMS              *ip,
                              char                     *szTitle,
                              PINChI2                  *pINChI2[INCHI_NUM],
                              PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                              INCHI_IOSTREAM           *inp_file,
                              INCHI_IOSTREAM           *log_file,
                              INCHI_IOSTREAM           *out_file,
                              INCHI_IOSTREAM           *prb_file,
                              ORIG_ATOM_DATA           *orig_inp_data,
                              ORIG_ATOM_DATA           *prep_inp_data,
                              long                     num_inp,
                              INCHI_IOS_STRING         *strbuf,
                              unsigned char            save_opt_bits,
                              char					 **sinchi,
                              char					 **saux);
int OAD_ProcessOneStructureNoEdits( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    int                      *n_pzz,
                                    char					 **sinchi,
                                    char					 **saux);
int OAD_ProcessOneStructure105Plus( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    char					 **sinchi,
                                    char					 **saux);
int ValidateAndPreparePolymerAndPseudoatoms(struct tagINCHI_CLOCK *ic,
                                            struct tagCANON_GLOBALS  *CG,
                                            STRUCT_DATA *sd,
                                            INPUT_PARMS *ip,
                                            char *szTitle,
                                            PINChI2 *pINChI2[INCHI_NUM],
                                            PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                                            INCHI_IOSTREAM *inp_file,
                                            INCHI_IOSTREAM *log_file,
                                            INCHI_IOSTREAM *out_file,
                                            INCHI_IOSTREAM *prb_file,
                                            ORIG_ATOM_DATA *orig_inp_data,
                                            ORIG_ATOM_DATA *prep_inp_data,
                                            long num_inp,
                                            INCHI_IOS_STRING *strbuf,
                                            unsigned char save_opt_bits,
                                            int *mind_polymers);
int CreateOneStructureINChI( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic,
                             STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                             PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                             int iINChI,
                             INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                             INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file,
                             ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                             COMP_ATOM_DATA composite_norm_data2[][TAUT_NUM + 1],
                             long num_inp, INCHI_IOS_STRING *strbuf, NORM_CANON_FLAGS *pncFlags );
int PreprocessOneStructure( struct tagINCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip,
                            ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data );
int RepeatedlyRenumberAtomsAndRecalcINCHI( struct tagINCHI_CLOCK *ic, CANON_GLOBALS *CG,
                                           STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                                           PINChI2 *pINChI[INCHI_NUM], PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                                           INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *plog,
                                           INCHI_IOSTREAM *pout, INCHI_IOSTREAM *pprb,
                                           ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                                           long *num_inp, STRUCT_FPTRS *pStructPtrs,
                                           int *nRet, int have_err_in_GetOneStructure,
                                           long *num_err, int output_error_inchi, INCHI_IOS_STRING *strbuf,
                                           unsigned long *pulTotalProcessingTime, char *pLF, char *pTAB,
                                           long int nrepeat);
int bIsStructChiral( PINChI2 *pINChI2[INCHI_NUM], int num_components[] );


/* ORIG_ATOM_DATA  */

int  OrigAtData_Duplicate( ORIG_ATOM_DATA *new_orig_atom, ORIG_ATOM_DATA *orig_atom );

int  OrigAtData_RemoveAtom(ORIG_ATOM_DATA *orig_at_data, int iatom);
int  OrigAtData_AddSingleStereolessBond( int this_atom, int other_atom,
                                         inp_ATOM *at, int *num_inp_bonds );
int  OrigAtData_AddBond( int this_atom, int other_atom, inp_ATOM *at,
                         int bond_type, int bond_stereo, int *num_bonds );
int  OrigAtData_RemoveBond( int this_atom, int other_atom, inp_ATOM *at,
                            int *bond_type, int *bond_stereo, int *num_inp_bonds );
int  OrigAtData_RemoveHalfBond( int this_atom, int other_atom, inp_ATOM *at,
                                int *bond_type, int *bond_stereo );
int  OrigAtData_IncreaseBondOrder( int this_atom, int other_atom, inp_ATOM *at );
int  OrigAtData_DecreaseBondOrder( int this_atom, int other_atom, inp_ATOM *at );
int  OrigAtData_SaveMolfile( ORIG_ATOM_DATA *orig_inp_data, STRUCT_DATA *sd,
                             INPUT_PARMS *ip, long num_inp, INCHI_IOSTREAM *out_file );
void OrigAtData_DebugTrace( ORIG_ATOM_DATA *at_data );
int OAD_StructureEdits_Apply( STRUCT_DATA *sd,
                             INPUT_PARMS *ip,
                             ORIG_ATOM_DATA *orig_at_data,
                             OAD_StructureEdits *ed,
                             int *ret);
void OAD_CollectFragmentBondsAndAtoms( ORIG_ATOM_DATA  *at_data,
                                       int nforbidden,		
                                       int *forbidden_orig,
                                       int *n_fragbonds,
                                       int **fragbonds,
                                       int *n_fragatoms,
                                       int *fragatoms,
                                       int *err,
                                       char *pStrErr);

void winchi_calc_inchikey(  int ret,
                            int            *ikflag,
                            INPUT_PARMS    *ip,
                            INCHI_IOSTREAM *out_file,
                            INCHI_IOSTREAM *log_file);

/* inp_ATOM */
int Inp_Atom_GetBondType(inp_ATOM *at, int iatom1, int iatom2);


/* ORIG_STRUCT */

int OrigStruct_FillOut( struct tagCANON_GLOBALS *pCG,
                        ORIG_ATOM_DATA *orig_inp_data,
                        ORIG_STRUCT *pOrigStruct,
                        STRUCT_DATA *sd );
void OrigStruct_Free( ORIG_STRUCT *pOrigStruct );
int     ReadWriteInChI( struct tagINCHI_CLOCK *ic,
                    struct tagCANON_GLOBALS *pCG,
                    INCHI_IOSTREAM *pInp,
                    INCHI_IOSTREAM *pOut,
                    INCHI_IOSTREAM *pLog,
                    INPUT_PARMS *ip_inp,
                    STRUCT_DATA *sd_inp,
                    /* the following are InChI library-specific parameters */
                    inp_ATOM **at,
                    int *num_at,
                    int *num_bonds,
                    OAD_Polymer **polymer,
                    OAD_V3000 **v3000,
                    /* end of InChI library-specific parameters */
                    char *szMsg,
                    int nMsgLen,
                    unsigned long WarningFlags[2][2] );
int CompareHillFormulasNoH( const char *f1, const char *f2,
                            int *num_H1, int *num_H2 );
int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data,
                             INP_ATOM_DATA2 *all_inp_norm_data,
                             int num_components );


/* POLYMERS */

void EditINCHI_HidePolymerZz( INCHI_IOSTREAM *out, int n_pzz, int n_zy );




/* MISCELLANEOUS */


void SplitTime( unsigned long ulTotalTime, int *hours, int *minutes, int *seconds, int *mseconds );
void set_line_separators( int bINChIOutputOptions, char **pLF, char **pTAB );
void save_command_line( int argc, char *argv[], INCHI_IOSTREAM *plog );
void emit_empty_inchi( INPUT_PARMS *ip, long num_inp,
                       char *pLF, char *pTAB, INCHI_IOSTREAM *pout );

#ifndef COMPILE_ANSI_ONLY
void eat_keyboard_input( void );
int user_quit( struct tagINCHI_CLOCK *ic, const char *msg, unsigned long ulMaxTime );
#endif



#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif    /* _ICHIMAIN_H_ */
