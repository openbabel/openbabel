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


#ifndef __INCHIMAIN_H__
#define __INCHIMAIN_H__

#define ESC_KEY       27
#define INCHI_SEGM_BUFLEN  64000

/********************************************************************/
typedef struct tagStructData {
    unsigned long ulStructTime;
    int           nErrorCode;
    int           nErrorType;
    int           nStructReadError;
    char          pStrErrStruct[STR_ERR_LEN];
    long          fPtrStart;  /* or number of processed structures */
    long          fPtrEnd;    /* or number of errors */
    int           bXmlStructStarted;
    int           bUserQuit;
    int           bUserQuitComponent;
    int           bUserQuitComponentDisplay;
    int           bChiralFlag;
    /* information related to normal or disconnected layers */
    int           num_taut[INCHI_NUM];
    int           num_non_taut[INCHI_NUM];
    INCHI_MODE     bTautFlags[INCHI_NUM];     /* reconnected does not have TG_FLAG_DISCONNECT_COORD_DONE flag */
    INCHI_MODE     bTautFlagsDone[INCHI_NUM]; /* reconnected does not have TG_FLAG_DISCONNECT_COORD_DONE flag */
    int           num_components[INCHI_NUM]; /* number of allocated INChI, INChI_Aux data structures */
    /* debugging info */
#if( bRELEASE_VERSION == 0 )
    int           bExtract;
#endif

} STRUCT_DATA;

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int ReadCommandLineParms( int argc, const char *argv[], INPUT_PARMS *ip, char *szSdfDataValue,
                          unsigned long *ulDisplTime, int bReleaseVersion, INCHI_FILE *log_file );
void HelpCommandLineParms( INCHI_FILE *f );
int OpenFiles( FILE **inp_file, FILE **output_file, FILE **log_file, FILE **prb_file, INPUT_PARMS *ip );
#ifndef INCHI_ANSI_ONLY
int DisplayStructure( inp_ATOM *at, int num_at, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H nNumRemovedProtonsIsotopic[],
                      int bIsotopic, int j /*bTautomeric*/, INChI **cur_INChI, INChI_Aux **cur_INChI_Aux,
                      int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle );
void FillTableParms( SET_DRAW_PARMS *sdp, INChI **cur_INChI, INChI_Aux **cur_INChI_Aux, INCHI_MODE nMode, int bShowIsotopic, int bShowTaut );
void FillCompositeTableParms( SET_DRAW_PARMS *sdp, AT_NUMB StereoFlags,
                     INCHI_MODE nMode, int bShowIsotopic, int bShowTaut );
#endif
int PrintInputParms( INCHI_FILE *log_file, INPUT_PARMS *ip );
const char *ErrMsg( int nErrorCode );
int SortAndPrintINChI( INCHI_FILE *output_file, char *pStr, int nStrLen, INCHI_FILE *log_file,
                      INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                      COMP_ATOM_DATA composite_norm_data[INCHI_NUM][TAUT_NUM+1],
                      ORIG_STRUCT *pOrigStruct, int num_components[INCHI_NUM],
                      int num_non_taut[INCHI_NUM], int num_taut[INCHI_NUM],
                      INCHI_MODE bTautFlags[INCHI_NUM], INCHI_MODE bTautFlagsDone[INCHI_NUM],
                      NORM_CANON_FLAGS *pncFlags, long num_inp,
                      PINChI2 *pINChI[INCHI_NUM], PINChI_Aux2 *pINChI_Aux[INCHI_NUM], int *bSortPrintINChIFlags );
void FreeAllINChIArrays( PINChI2 *pINChI[INCHI_NUM], PINChI_Aux2 *pINChI_Aux[INCHI_NUM], int num_components[2] );
void FreeINChIArrays( PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components );
void SplitTime( unsigned long ulTotalTime, int *hours, int *minutes, int *seconds, int *mseconds );

int ReadTheStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, FILE *inp_file, ORIG_ATOM_DATA *orig_inp_data,
                     int inp_index, int *out_index );
int TreatReadTheStructureErrors(  STRUCT_DATA *sd, INPUT_PARMS *ip, int nLogMask, 
                                  FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                                  FILE *prb_file, /*^^^ was: INCHI_FILE */
                                  ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen );

int GetOneComponent( STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_FILE *log_file, INCHI_FILE *output_file,
                     INP_ATOM_DATA *inp_cur_data,
                     ORIG_ATOM_DATA *orig_inp_data, int i, long num_inp, char *pStr, int nStrLen );
int CreateOneComponentINChI( STRUCT_DATA *sd, INPUT_PARMS *ip, INP_ATOM_DATA *inp_cur_data, ORIG_ATOM_DATA *orig_inp_data,
                            PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int iINChI,
                            int i, long num_inp, INP_ATOM_DATA **inp_norm_data,
                            NORM_CANON_FLAGS *pncFlags, INCHI_FILE *log_file );
int TreatCreateOneComponentINChIError(STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data,
                                     int i, long num_inp,
                                     FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                                     FILE *prb_file, /*^^^ was: INCHI_FILE */
                                     char *pStr, int nStrLen );
int TreatCreateINChIWarning(STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, long num_inp,
                                     FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                                     FILE *prb_file, /*^^^ was: INCHI_FILE */
                                     char *pStr, int nStrLen );

#if( TEST_RENUMB_ATOMS == 1 || READ_INCHI_STRING == 1 ) /*  { */
int CompareINChI( INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2 );
#endif

void eat_keyboard_input( void );
int user_quit( const char *msg, unsigned long ulMaxTime );

int GetOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                     FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                     FILE *prb_file, /*^^^ was: INCHI_FILE */
                     ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen, STRUCT_FPTRS *struct_fptrs );
int ProcessOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                         PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                         FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                         FILE *prb_file, /*^^^ was: INCHI_FILE */
                         ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                         long num_inp, char *pStr, int nStrLen ); 

int CreateOneStructureINChI( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                         PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM], int iINChI,
                         FILE *inp_file, INCHI_FILE *log_file, INCHI_FILE *output_file, 
                         FILE *prb_file, /*^^^ was: INCHI_FILE */
                         ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                         COMP_ATOM_DATA composite_norm_data2[][TAUT_NUM+1],
                         long num_inp, char *pStr, int nStrLen, NORM_CANON_FLAGS *pncFlags );

int bIsStructChiral( PINChI2 *pINChI2[INCHI_NUM], int num_components[] );

int PreprocessOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data );

int FillOutOrigStruct( ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, STRUCT_DATA *sd );
void FreeOrigStruct(  ORIG_STRUCT *pOrigStruct);


int ReadWriteInChI(
#if ( defined(BUILD_CINCHI_WITH_INCHIKEY) && !defined(INCHI_LIB) )                   
                   FILE *pInp, FILE *output_file, FILE *log_file,
#else
                   INCHI_FILE *pInp, 
#endif
                   INCHI_FILE *pOut, INCHI_FILE *pLog,
                   INPUT_PARMS *ip,  STRUCT_DATA *sd, inp_ATOM **at, int *num_at,
                   char *szMsg, int nMsgLen, unsigned long WarningFlags[2][2] );

int CompareHillFormulasNoH( const char *f1, const char *f2, int *num_H1, int *num_H2 );


#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /*  __INCHIMAIN_H__ */
