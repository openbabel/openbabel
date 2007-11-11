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


#ifndef __INCHIDRP_H__
#define __INCHIDRP_H__

#ifndef INCHI_ANSI_ONLY /* { */
/********************************************
 * Parameters for the structure drawing
 ********************************************/
#define TDP_LEN_LBL      16  /* length of a label (label: Req., Shown, Found) */
/* #define TDP_NUM_LBL 3  */ /* number of labels */
/* #define TDP_NUM_PAR 3  */ /* number of types per label (types: B/T, I/N, S) */
typedef enum tagTblTypes {itBASIC, itISOTOPIC, itSTEREO, TDP_NUM_PAR} TBL_TYPES; /*  types */
typedef enum tagTblLabels{ ilSHOWN,  TDP_NUM_LBL} TBL_LABELS; /*  labels */
typedef struct tagTblDrawPatms {
    char   ReqShownFoundTxt[TDP_NUM_LBL][TDP_LEN_LBL];
    char   ReqShownFound[TDP_NUM_LBL][TDP_NUM_PAR];
    int    nOrientation;  /* 10*degrees: 0 or 2700 */
    int    bDrawTbl;
} TBL_DRAW_PARMS;
/*********************************************/
typedef struct tagDrawParmsSettings {
    TBL_DRAW_PARMS *tdp;
    unsigned long  ulDisplTime;
    int            bOrigAtom;
    int            nFontSize;
} SET_DRAW_PARMS;  /* input only: how to draw or calculate */
/*********************************************/
typedef struct tagReturnedDrawParms {
    int       bEsc;
} RET_DRAW_PARMS;
/*********************************************/
typedef struct tagPersistDrawParms {
    int rcPict[4];
} PER_DRAW_PARMS; /* saved between displaying different structures */
/*********************************************/
typedef struct tagDrawParms {
    SET_DRAW_PARMS  sdp;   /* how to draw: fill on the 1st call */
    RET_DRAW_PARMS  rdp;   /* returned when drawing window is closed */
    PER_DRAW_PARMS *pdp;   /* persistent: save between calls (window size) */
#ifndef INCHI_LIB
#ifndef INCHI_ANSI_ONLY
    AT_NUMB   *nEquLabels; /* num_inp_atoms elements, value>0 marks atoms in the set #value  */
    AT_NUMB    nNumEquSets;  /* max mark value */
    AT_NUMB    nCurEquLabel; /* current mark */
#endif
#endif
} DRAW_PARMS; /* Settings: How to draw the structure */

#endif /* } INCHI_ANSI_ONLY */

#if( ACD_LABS_VERSION == 1 )
#define MAX_NUM_PATHS 5
#else
#define MAX_NUM_PATHS 4
#endif

typedef enum tagInputType { INPUT_NONE=0, INPUT_MOLFILE=1, INPUT_SDFILE=2, INPUT_INCHI_XML=3, INPUT_INCHI_PLAIN=4, INPUT_CMLFILE=5, INPUT_INCHI=6, INPUT_MAX } INPUT_TYPE;
typedef struct tagInputParms {
    char            szSdfDataHeader[MAX_SDF_HEADER+1];
    char           *pSdfLabel;
    char           *pSdfValue;
    long            lSdfId;
    long            lMolfileNumber;
#ifndef INCHI_ANSI_ONLY
    DRAW_PARMS      dp;
    PER_DRAW_PARMS  pdp;
    TBL_DRAW_PARMS  tdp;
#endif
/*
  -- Files --
  ip->path[0] => Input
  ip->path[1] => Output (INChI)
  ip->path[2] => Log
  ip->path[3] => Problem structures
  ip->path[4] => Errors file (ACD)

*/
    const char     *path[MAX_NUM_PATHS];
    int             num_paths;
    long            first_struct_number;
    long            last_struct_number;
    INPUT_TYPE      nInputType;
    INCHI_MODE      nMode;
    int             bAbcNumbers;
    /*int             bXml;*/
    int             bINChIOutputOptions; /* !(ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) */
    int             bCtPredecessors;
    int             bXmlStarted;
    int             bDisplayEachComponentINChI;

    long            msec_MaxTime;   /* was ulMaxTime; max time to run ProsessOneStructure */
    long            msec_LeftTime;

    long            ulDisplTime; /* not used: max structure or question display time */
    int             bDisplay;
    int             bDisplayIfRestoreWarnings; /* InChI->Struct debug */
    int             bMergeAllInputStructures;
    int             bSaveWarningStructsAsProblem;
    int             bSaveAllGoodStructsAsProblem;
    int             bGetSdfileId;
    int             bGetMolfileNumber;  /* read molfile number from the name line like "Structure #22" */
    int             bCompareComponents; /* see flags CMP_COMPONENTS, etc. */
    int             bDisplayCompositeResults;
    int             bDoNotAddH;
    int             bNoStructLabels;
    int             bChiralFlag;
    int             bAllowEmptyStructure;
	/*^^^ */
    int             bCalcInChIKey;
	/*^^^ */
    INCHI_MODE      bTautFlags;
    INCHI_MODE      bTautFlagsDone;

#if( READ_INCHI_STRING == 1 )
    int             bReadInChIOptions;
#endif

/* post v.1 features */
#if( UNDERIVATIZE == 1 )
    int             bUnderivatize;
#endif
#if( RING2CHAIN == 1 )
    int             bRing2Chain;
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    int             bIngnoreUnchanged;
#endif

} INPUT_PARMS;

#endif /* __INCHIDRP_H__ */
