/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* #define CHECK_WIN32_VC_HEAP */
#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichicomp.h"
#include "ichi.h"
#include "ichitime.h"
#include "util.h"
#include "strutil.h"
#include "ichi_io.h"

/* reverse InChI */
#include "ichimain.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichister.h"
#include "strutil.h"
#include "ichisize.h"
#include "ichiring.h"
#include "ichinorm.h"
#include "ichierr.h"

#include "ichirvrs.h"


/*^^^ */
#if ( defined(TARGET_API_LIB) || defined(TARGET_EXE_STANDALONE) )
#include "inchi_api.h"
#endif
/*^^^ */



typedef struct tagLine 
{
    char *str;
    int   len;
    int   len_alloc;
    int   c;
} SEGM_LINE;
#define SEGM_LINE_ADD 128

typedef struct tagOneLinkedBond 
{
    AT_NUMB neigh; /* canonical number of a neighbor */
    AT_NUMB prev;  /* position of the previous neighbor in the list */
} ONE_LINKED_BOND;

typedef struct tagLinkedBonds 
{
    ONE_LINKED_BOND *pBond;
    int len;
    int len_alloc;
}LINKED_BONDS;
#define LINKED_BOND_ADD  128

typedef enum tagModeProtonIsoExchgH 
{
    MODE_PIXH_UNDEFINED,     /* 0 */
    MODE_PIXH_ADD_TO_FIRST,  /* 1 */
    MODE_PIXH_ADD_TO_EACH,   /* 2 */
    MODE_PIXH_ADD_A_PIXH_COMPONENT, /* 3 */
    MODE_PIXH_KEEP_TOTALS         /* 4 */
} MODE_PIXH;


/* local prototypes */
static int GetInChIFormulaNumH(INChI *pInChI, int *nNumH);
static int GetInChINumH(INChI *pInChI, int *nNumH);
static int GetInChIIsoH(INChI *pInChI, int nNumIsotopicH[NUM_H_ISOTOPES]);

static int getInChIChar(INCHI_IOSTREAM *pInp);
static int AddInChIChar(INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken);
static int AddLinkedBond(AT_NUMB at1, AT_NUMB at2, AT_NUMB num_at, LINKED_BONDS *pLB);
static int bInChIHasReconnectedMetal(INChI *pInChI);
static int SetProtonsAndXchgIsoH(int bInChI2Structure, 
                                 int bReqSplitOutputInChI, 
                                 int bReqProtonsForEachComponent,
                                 int bReqNonTaut, int bReqStereo, 
                                 int num_components[INCHI_NUM],
                                 MODE_PIXH nModeProtonIsoExchgH[INCHI_NUM], 
                                 InpInChI *OneInput);
#if ( FIX_DALKE_BUGS == 1 )
static int SetHillFormFromInChI(InpInChI *OneInput);
#endif

static int nGetInChISegment(INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken);

static int CopySegment(INChI *pInChITo, INChI *pInChIFrom, int StereoType, 
                       int bIsotopicTo, int bIsotopicFrom);
static int nFillOutProtonMobileH(INChI *pInChI);
static int nProtonCopyIsotopicInfo(INChI *pInChI_to, INChI *pInChI_from);
static int CopyAtomNumbers(INChI *pInChI_To, int bIsoTo, INChI *pInChI_From, int bIsoFrom);

static int ParseSegmentFormula(const char *str, int bMobileH, INChI *pInpInChI[], 
                               int nNumComponents[]);
static int ParseSegmentConnections(const char *str, int bMobileH, INChI **pInpInChI, 
                                   int *pnNumComponents, int *pbAbc);
static int ParseSegmentMobileH(const char *str, int bMobileH, 
                               INChI *pInpInChI[], int pnNumComponents[], int *pbAbc);
static int ParseSegmentCharge(const char *str, int bMobileH, 
                              INChI *pInpInChI[], int nNumComponents[]);
static int ParseSegmentProtons(const char *str, int bMobileH, 
                               REM_PROTONS nNumProtons[], int nNumComponents[]);
static int ParseSegmentSp2(const char *str, int bMobileH, 
                           INChI *pInpInChI[], int nNumComponents[], int state, int *pbAbc);
static int ParseSegmentSp3(const char *str, int bMobileH, 
                           INChI *pInpInChI[], int nNumComponents[], int state, int *pbAbc);
static int ParseSegmentSp3m(const char *str, int bMobileH, 
                            INChI *pInpInChI[], int nNumComponents[], int state);
static int bIsSp3LayerNotEmpty(INChI *pInpInChI[], int bMobileH, 
                               int bIso, int nNumComponents);
static int ParseSegmentSp3s(const char *str, int bMobileH, 
                            INChI *pInpInChI[], int s[TAUT_NUM][2], int ppnNumComponents[], int state);
static int ParseSegmentIsoAtoms(const char *str, int bMobileH, INChI *pInpInChI[], 
                                int nNumComponents[], int state, int *pbAbc);
static int ParseSegmentIsoExchgH(const char *str, int bMobileH, REM_PROTONS nNumProtons[], 
                                 int nNumComponents[], int state, int *pbAbc);
static int ParseSegmentPerm(const char *str, int bMobileH, INChI *pInpInChI[], 
                            int ppnNumComponents[], int state, int *pbAbc);
#if ( FIX_ISO_FIXEDH_BUG_READ == 1 )
static int bIsoMayBeArranged(int bInchi2Struct, int iso_diff[NUM_H_ISOTOPES], 
                             REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                             INChI *pInpInChI[INCHI_NUM][TAUT_NUM], int nNumComponents[INCHI_NUM][TAUT_NUM], int iINChI);
#endif

static int ReadInChILine(INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, 
                         char **pStr, int *pState,
                         INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                         int nNumComponents[INCHI_NUM][TAUT_NUM],
                         REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                         int s[INCHI_NUM][TAUT_NUM][2],
                         int *bStdFormat, 
                         int *bInputHasSaveOpt, unsigned char *inp_save_opt_bits);

int InChILine2Data(INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, 
                   char **pStr, int *pState, int *nErr,
                   INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                   int nNumComponents[INCHI_NUM][TAUT_NUM],
                   REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                   int s[INCHI_NUM][TAUT_NUM][2], 
                   int bReadCoord, int bInchi2Struct, 
                   INCHI_MODE nMode,
                   int *bStdFormat, 
                   int *bInputHasSaveOpt, unsigned char *inp_save_opt_bits);

static int ReadInChICoord(INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, int *pState,
                          INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                          int nNumComponents[INCHI_NUM][TAUT_NUM]);

static int OutputInChIAsRequested(INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, 
                                  ICHICONST INPUT_PARMS *ip_inp, 
                                  STRUCT_DATA *sd_inp, 
                                  InpInChI *OneInput, 
                                  int num_components[INCHI_NUM],
                                  MODE_PIXH nModeProtonIsoExchgH[INCHI_NUM], 
                                  long num_inp, unsigned char save_opt_bits);

static int ParseAuxSegmentVersion(const char *str, int bMobileH, INChI *pInpInChI[], 
                                  int ppnNumComponents[], int state);
static int ParseAuxSegmentNumbers(const char *str, int bMobileH, INChI *pInpInChI[], 
                                  int ppnNumComponents[], int state, int *pbAbc);
static int ParseAuxSegmentAtomEqu(const char *str, int bMobileH, INChI *pInpInChI[], 
                                  int ppnNumComponents[], int state);
static int ParseAuxSegmentGroupEqu(const char *str, int bMobileH, INChI *pInpInChI[], 
                                   int ppnNumComponents[], int state);
static int ParseAuxSegmentSp3Inv(const char *str, int bMobileH, INChI *pInpInChI[], 
                                 int ppnNumComponents[], int state);
static int ParseAuxSegmentSp3InvNumbers(const char *str, int bMobileH, INChI *pInpInChI[], 
                                        int ppnNumComponents[], int state);
static int ParseAuxSegmentReverseCRV(const char *str, int bMobileH, INChI *pInpInChI[], 
                                     int ppnNumComponents[], int state);
static int ParseAuxSegmentReverseAtoms(const char *str, int bMobileH, INChI *pInpInChI[], 
                                       int ppnNumComponents[], int state);
static int ParseAuxSegmentReverseBonds(const char *str, int bMobileH, INChI *pInpInChI[], 
                                       int ppnNumComponents[], int state);
static int ParseAuxSegmentReverseXYZ(const char *str, int bMobileH, XYZ_COORD **ppXYZ, 
                                     INChI *pInpInChI[], int ppnNumComponents[], int state);
static int AddAuxSegmentCoord(int nRet, XYZ_COORD *pXYZ, int nLenXYZ, 
                              INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                              int nNumComponents[INCHI_NUM][TAUT_NUM]);

static const char *getInchiStateReadErr(int stat);
static const char *getInchiErrName(int nErr);

#define SEG_END '/'
/* the following 2 definitions are used to allow tab-delimited InChI input - 2008-11-17 DT */
#define INCHI_INP_EOL(X) ((X)=='\n' || (X)=='\r' || (X)=='\t')
/*#define INCHI_TOKEN "/\n\r\t"*/
#define INCHI_TOKEN "/\n\r\t\\"

typedef enum tagInChI_STATE 
{
    /* M */
    IST_MOBILE_H_FORMULA,          /* 0 */
    IST_MOBILE_H_CONNECTIONS,      /* 1 */
    IST_MOBILE_H,                  /* 2 */
    IST_MOBILE_H_CHARGE,           /* 3 */
    IST_MOBILE_H_PROTONS,          /* 4 */
    IST_MOBILE_H_SP2,              /* 5 */
    IST_MOBILE_H_SP3,              /* 6 */
    IST_MOBILE_H_SP3_M,            /* 7 */
    IST_MOBILE_H_SP3_S,            /* 8 */

    /* Fork */
    IST_MOBILE_H_ISO_LAYER_FORK,   /* 9 */

    /* MI */
    IST_MOBILE_H_ISO_ATOMS,        /* 10 */
    IST_MOBILE_H_ISO_EXCH_H,       /* 11 */
    IST_MOBILE_H_ISO_SP2,          /* 12 */
    IST_MOBILE_H_ISO_SP3,          /* 13 */
    IST_MOBILE_H_ISO_SP3_M,        /* 14 */
    IST_MOBILE_H_ISO_SP3_S,        /* 15 */

    /* Fork */
    IST_FIXED_H_LAYER_FORK,        /* 16 */
    
    /* F */
    IST_FIXED_H_FORMULA,           /* 17 */
    IST_FIXED_H,                   /* 18 */
    IST_FIXED_H_CHARGE,            /* 19 */
    IST_FIXED_H_SP2,               /* 20 */
    IST_FIXED_H_SP3,               /* 21 */
    IST_FIXED_H_SP3_M,             /* 22 */
    IST_FIXED_H_SP3_S,             /* 23 */
    IST_FIXED_H_PERMUTATION,       /* 24 */

    /* Fork */
    IST_FIXED_H_ISO_LAYER_FORK,    /* 25 */

    /* FI */
    IST_FIXED_H_ISO_ATOMS,         /* 26 */
    IST_FIXED_H_ISO_LAYER,         /* 27 */
    IST_FIXED_H_ISO_SP2,           /* 28 */
    IST_FIXED_H_ISO_SP3,           /* 29 */
    IST_FIXED_H_ISO_SP3_M,         /* 30 */
    IST_FIXED_H_ISO_SP3_S,         /* 31 */
    IST_FIXED_H_ISO_PERMUTATION,   /* 32 */

    /* Reconnected */
    IST_RECONNECTED_LAYER_FORK,    /* 33 */
    IST_RECONNECTED_FORMULA,       /* 34 */

    /* Other reading errors */
    IST_MATERIAL_BALANCE_ERROR,    /* 35 */

    IST_END                     =    -1
}INCHI_STATE;

#define IST_HAPPENED_IN_RECMET   100

typedef struct tagInchiReadErrMsg 
{
    int         stat;
    const char  *msg;
} INCHI_READ_ERR_MSG;

ICHICONST INCHI_READ_ERR_MSG irErrMsg[] = 
{
    /* M */
   {IST_MOBILE_H_FORMULA,            "MOBILE_H_FORMULA"          },
   {IST_MOBILE_H_CONNECTIONS,        "MOBILE_H_CONNECTIONS"      },
   {IST_MOBILE_H,                    "MOBILE_H"                  },
   {IST_MOBILE_H_CHARGE,             "MOBILE_H_CHARGE"           },
   {IST_MOBILE_H_PROTONS,            "MOBILE_H_PROTONS"          },
   {IST_MOBILE_H_SP2,                "MOBILE_H_SP2"              },
   {IST_MOBILE_H_SP3,                "MOBILE_H_SP3"              },
   {IST_MOBILE_H_SP3_M,              "MOBILE_H_SP3_/m"           },
   {IST_MOBILE_H_SP3_S,              "MOBILE_H_SP3_/s"           },

    /* Fork */
   {IST_MOBILE_H_ISO_LAYER_FORK,     "MOBILE_H_ISO_LAYER_FORK"   },

    /* MI */
   {IST_MOBILE_H_ISO_ATOMS,          "MOBILE_H_ISO_ATOMS"        },
   {IST_MOBILE_H_ISO_EXCH_H,         "MOBILE_H_ISO_EXCH_H"       },
   {IST_MOBILE_H_ISO_SP2,            "MOBILE_H_ISO_SP2"          },
   {IST_MOBILE_H_ISO_SP3,            "MOBILE_H_ISO_SP3"          },
   {IST_MOBILE_H_ISO_SP3_M,          "MOBILE_H_ISO_SP3_/m"       },
   {IST_MOBILE_H_ISO_SP3_S,          "MOBILE_H_ISO_SP3_/s"       },

    /* Fork */
  {IST_FIXED_H_LAYER_FORK,           "FIXED_H_LAYER_FORK"        },

    /* F */
   {IST_FIXED_H_FORMULA,             "FIXED_H_FORMULA"           },
   {IST_FIXED_H,                     "FIXED_H"                   },
   {IST_FIXED_H_CHARGE,              "FIXED_H_CHARGE"            },
   {IST_FIXED_H_SP2,                 "FIXED_H_SP2"               },
   {IST_FIXED_H_SP3,                 "FIXED_H_SP3"               },
   {IST_FIXED_H_SP3_M,               "FIXED_H_SP3_/m"            },
   {IST_FIXED_H_SP3_S,               "FIXED_H_SP3_/s"            },
   {IST_FIXED_H_PERMUTATION,         "FIXED_H_PERMUTATION"       },

    /* Fork */
   {IST_FIXED_H_ISO_LAYER_FORK,      "FIXED_H_ISO_LAYER_FORK"    },

    /* FI */
   {IST_FIXED_H_ISO_ATOMS,           "FIXED_H_ISO_ATOMS"         },
   {IST_FIXED_H_ISO_LAYER,           "FIXED_H_ISO_LAYER"         },
   {IST_FIXED_H_ISO_SP2,             "FIXED_H_ISO_SP2"           },
   {IST_FIXED_H_ISO_SP3,             "FIXED_H_ISO_SP3"           },
   {IST_FIXED_H_ISO_SP3_M,           "FIXED_H_ISO_SP3_m"         },
   {IST_FIXED_H_ISO_SP3_S,           "FIXED_H_ISO_SP3_s"         },
   {IST_FIXED_H_ISO_PERMUTATION,     "FIXED_H_ISO_PERMUTATION"   },

    /* Reconnected */
   {IST_RECONNECTED_LAYER_FORK,      "RECONNECTED_LAYER_FORK"    },
   {IST_RECONNECTED_FORMULA,         "RECONNECTED_FORMULA"       },

   {IST_MATERIAL_BALANCE_ERROR,      "MATERIAL_BALANCE"          },

   {IST_END,                         "Unknown Error"             }
};



typedef enum tagCopySegmentType 
{
    CPY_SP2,
    CPY_SP3,
    CPY_SP3_M,
    CPY_SP3_S,
    CPY_ISO_AT
} COPY_SEG_TYPE;

#define NSTRLEN 64000
#define MAX_MSG_LEN 512
#define MAX_MSG_BUF_LEN 128


/*************************************************************************************/
const char *getInchiStateReadErr(int stat)
{
    int i, bRecMet = 0;
    static char szMsg[128];
    if ( stat >= IST_HAPPENED_IN_RECMET ) 
    {
        bRecMet = 1;
        stat -= IST_HAPPENED_IN_RECMET;
    }
    for ( i = 0; 0 <= irErrMsg[i].stat && stat != irErrMsg[i].stat; i ++ )
        ;
    sprintf(szMsg, 
#if ( FIX_DALKE_BUGS == 1 )
        "%s%.100s",
#else
        "%s%s",
#endif
        irErrMsg[i].msg, bRecMet? ", Reconnected layer" : "");
    return szMsg;
}

/**************************************************************************************/
const char *getInchiErrName(int nErr)
{
    switch ( nErr ) 
    {
    case RI_ERR_ALLOC:
        return "Allocation failed";
    case RI_ERR_PROGR:
        return "Program error";
    case RI_ERR_SYNTAX:
        return "Syntax error";
    case RI_ERR_EOL:
        return "End of line";
    }
    return "Unknown error";
}
#if ( FIX_DALKE_BUGS == 1 )
/*****************************************************************************************/
int SetHillFormFromInChI(InpInChI *OneInput)
{
    int iINChI, iTaut, iComp, num_diff;
    INChI *pINChI;
    char *szHillFormulaOld;
    for ( iINChI = 0, num_diff = 0; iINChI < INCHI_NUM; iINChI ++ ) 
    {
        for ( iTaut = TAUT_NON; iTaut < TAUT_NUM; iTaut ++ ) 
        {
            for ( iComp = 0; iComp < OneInput->nNumComponents[iINChI][iTaut]; iComp ++ ) 
            {
                pINChI = &OneInput->pInpInChI[iINChI][iTaut][iComp];
                if ( !pINChI->nNumberOfAtoms || pINChI->bDeleted || !pINChI->szHillFormula || !pINChI->szHillFormula[0] ) 
                {
                    continue;
                }
                szHillFormulaOld = pINChI->szHillFormula;
                pINChI->szHillFormula = AllocateAndFillHillFormula(pINChI);
                num_diff += !pINChI->szHillFormula || !pINChI->szHillFormula[0] || strcmp(pINChI->szHillFormula, szHillFormulaOld);
                inchi_free(szHillFormulaOld);
            }
        }
    }
    return num_diff;
}
#endif



/********************** main entry point **********************************************/

int ReadWriteInChI(INCHI_IOSTREAM *pInp, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog,
                   INPUT_PARMS *ip_inp,  
                   STRUCT_DATA *sd_inp,
                   /* the following are InChI library-specific parameters */
                   inp_ATOM **at, int *num_at,
                   char *szMsg, int nMsgLen, unsigned long WarningFlags[2][2])
{                
    InpInChI OneInput;
    int      i, j, nReadStatus, ret, nErr, iINChI;
    char    *strHdr=NULL;
    char    *szCurHdr = NULL;
    int num_components[INCHI_NUM];
    int bReqNonTaut = (0 != ((ip_inp->nMode & REQ_MODE_BASIC) && 
                             (ip_inp->nMode & REQ_MODE_TAUT)));
    /*
    int bReqRecmet  = (0 != ((ip->bTautFlags & TG_FLAG_RECONNECT_COORD) &&
                             (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD)));
    */
    int bReqStereo  = (0 != (ip_inp->nMode & REQ_MODE_STEREO));
    int bHasSomeReconnected = 0, bHasSomeFixedH = 0, bHasMetal = 0;
    int nModeFlagsStereo = 0, bTautFlags = 0; /* InChI creation flags modifications derived from current InChI */
    MODE_PIXH nModeProtonIsoExchgH[INCHI_NUM];

    NORM_CANON_FLAGS ncFlags;
    NORM_CANON_FLAGS *pncFlags = &ncFlags;
    INPUT_PARMS ip_cur, *ip;
    STRUCT_DATA sd_cur, *sd;
    int  nMessageLen = MAX_MSG_LEN;
    char szMessage[MAX_MSG_LEN];
    int  nInitLenMessage;

    int  pState, bStereoType;
    int  bReqProtonsForEachComponent = 0;
    int  bReqSplitOutputInChI = 0;
    SEGM_LINE Line;
    SEGM_LINE *pLine = &Line;
    long          ulProcessingTime = 0;
    inchiTime     ulTStart;
    long          num_processed = 0, num_errors = 0;
    int  bPlainTabbedOutput;
    const char *pTAB;

#ifdef TARGET_API_LIB
    const int       bInChI2Structure = 0 != (ip_inp->bReadInChIOptions & READ_INCHI_TO_STRUCTURE);
    const int       bInChI2InChI     = 0 != (ip_inp->bReadInChIOptions & READ_INCHI_OUTPUT_INCHI);
#else
    const int       bInChI2Structure = 0 != (ip_inp->bReadInChIOptions & READ_INCHI_TO_STRUCTURE);
    const int       bInChI2InChI     = 0 != (ip_inp->bReadInChIOptions & READ_INCHI_OUTPUT_INCHI);
#endif
    const int       bReadCoord       = bInChI2Structure;
    long      num_inp=0;
   
    int bInputInStdFormat=0;
    int bInputHasSaveOpt=0;
    unsigned char inp_save_opt_bits=0;
    unsigned char save_opt_bits=0;

    ret         = 0;
    nReadStatus = RI_ERR_EOL;
    
    memset(szMessage, 0, sizeof(szMessage));
    memset(&OneInput, 0, sizeof(OneInput));
    memset(pLine, 0, sizeof(pLine[0]));
    if ( szMsg ) 
        szMsg[0] = '\0';


    while( nReadStatus != RI_ERR_EOF ) 
    {

        for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) 
        {
            for ( j = 0; j < TAUT_NUM; j ++ ) 
            {
                if ( OneInput.nNumProtons[iINChI][j].pNumProtons ) 
                {
                    inchi_free(OneInput.nNumProtons[iINChI][j].pNumProtons);
                    OneInput.nNumProtons[iINChI][j].pNumProtons = NULL;
                }
            }
        }
        
        memset(&OneInput, 0, sizeof(OneInput));
        memset(pncFlags, 0, sizeof(*pncFlags));        
        bStereoType = 0;
        ip_cur = *ip_inp;
        ip     = &ip_cur;
        sd_cur = *sd_inp;
        sd     = &sd_cur;
        bReqSplitOutputInChI          = 0 != (ip->bReadInChIOptions & READ_INCHI_SPLIT_OUTPUT);
        bReqProtonsForEachComponent   = bReqSplitOutputInChI && 
                                        0 != (READ_INCHI_KEEP_BALANCE_P & ip->bReadInChIOptions);

        bPlainTabbedOutput     = 0 != (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT);
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )    
        pTAB                   = bPlainTabbedOutput? "\t" : "\n";
#else
        pTAB                   = "\n";
#endif


        if ( bInChI2Structure ) 
        {
#if ( bRELEASE_VERSION == 1 )            
            bReqNonTaut = 1; /* bReqNonTaut=0 ignores Fixed-H layer in input InChI, for testing only */
#endif
            /* bReqRecmet  = 1; */
            bReqStereo  = 1;
            bReqSplitOutputInChI = 1;
            bReqProtonsForEachComponent = bReqNonTaut;
            ip->bTautFlags |= (TG_FLAG_DISCONNECT_COORD | TG_FLAG_RECONNECT_COORD);
            ip->nMode      |= (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_ISO);
        }
        
        
        /************************************************************************/
        /*                     Read InChI string                                */
        /************************************************************************/
        InchiTimeGet(&ulTStart);


        nReadStatus = InChILine2Data(pInp, pLine, &strHdr, &pState, &nErr, OneInput.pInpInChI,
                                    OneInput.nNumComponents, OneInput.nNumProtons,
                                    OneInput.s, bReadCoord, bInChI2Structure, 
                                    ip_inp->nMode, 
                                    &bInputInStdFormat,&bInputHasSaveOpt,&inp_save_opt_bits);


        
        ulProcessingTime += InchiTimeElapsed(&ulTStart);


        if ( (nReadStatus == RI_ERR_EOL || nReadStatus == RI_ERR_EOF) && !nErr &&
                               OneInput.nNumComponents[INCHI_BAS][TAUT_YES] 
                             + OneInput.nNumComponents[INCHI_BAS][TAUT_NON] ) 
        {
            /* InChI has been successfully read */
        
            ret = 0;
            num_inp ++;
            bHasSomeReconnected = 0;
            bHasSomeFixedH      = 0;

            /* Does not allow conversion non-standard->standard */
            /* (force target to be non-standard also)           */
            if ( ip_inp->bINChIOutputOptions & INCHI_OUT_STDINCHI ) 
            {
                if ( !bInputInStdFormat )   /* Input InChI is a non-standard one  */
                {
                    ip->bINChIOutputOptions &= ~INCHI_OUT_STDINCHI;
                    if ( szCurHdr && szCurHdr[0] ) 
                        inchi_ios_eprint( pLog, "Warning: forced conversion to non-standard InChI for non-std input, %s\n", szCurHdr );
                    else 
                        inchi_ios_eprint( pLog, "Warning: forced conversion to non-standard InChI for non-std input, Structure %ld\n", num_inp );
                }
            }            
            
            if ( ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT )
            {
                if (!bInputHasSaveOpt)
                {
                    /* Does not allow to create SaveOpt if the source lacks appendix */
                    ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
                    if ( szCurHdr && szCurHdr[0] ) 
                        inchi_ios_eprint( pLog, "Warning: ignore SaveOpt request for SaveOpt-less input, %s\n", szCurHdr );
                    else 
                        inchi_ios_eprint( pLog, "Warning: ignore SaveOpt request for SaveOpt-less input, Structure %ld\n", num_inp );
                }
                else
                {
                    /* Analyze existing and prepare new SaveOpt appendix */

                    if ( 0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD) )
                    {
                        /* RecMet requested */
                        if ( 0 != (inp_save_opt_bits & SAVE_OPT_RECMET) ) 
                        {
                            save_opt_bits |= SAVE_OPT_RECMET;
                        }
                        else
                        {
                            ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
                            if ( szCurHdr && szCurHdr[0] ) 
                                inchi_ios_eprint( pLog, "Warning: input created w/o RecMet - ignoring RecMet request, %s\n", szCurHdr );
                            else 
                                inchi_ios_eprint( pLog, "Warning: input created w/o RecMet - ignoring RecMet request, Structure %ld\n", num_inp );
                        }
                    }
                    if ( 0 != (ip->nMode & REQ_MODE_BASIC) )
                    {
                        /* FixedH requested */
                        if ( 0 != (inp_save_opt_bits & SAVE_OPT_FIXEDH) )
                            save_opt_bits |= SAVE_OPT_FIXEDH;
                        else
                        {
                            ip->nMode &= ~REQ_MODE_BASIC;
                            if ( szCurHdr && szCurHdr[0] ) 
                                inchi_ios_eprint( pLog, "Warning: input created w/o FixedH - ignoring FixedH request, %s\n", szCurHdr );
                            else 
                                inchi_ios_eprint( pLog, "Warning: input created w/o FixedH - ignoring FixedH request, Structure %ld\n", num_inp );
                        }
                    }
                    /* Copy from source SaveOpt the bits which we do not touch  */
                    /* while converting InChI:         SUU SLUUD KET 15T        */
                    if ( 0 != ( inp_save_opt_bits & SAVE_OPT_SUU) ) 
                        save_opt_bits |= SAVE_OPT_SUU;
                    if ( 0 != ( inp_save_opt_bits & SAVE_OPT_SLUUD) ) 
                        save_opt_bits |= SAVE_OPT_SLUUD;
                    if ( 0 != ( inp_save_opt_bits & SAVE_OPT_KET) ) 
                        save_opt_bits |= SAVE_OPT_KET;
                    if ( 0 != ( inp_save_opt_bits & SAVE_OPT_15T) ) 
                        save_opt_bits |= SAVE_OPT_15T;
                    /* Check if /SNon requested and turn OFF stereo bits if so */
                    if ( ! (ip->nMode & REQ_MODE_STEREO) )
                    {
                        save_opt_bits &= ~SAVE_OPT_SUU;
                        save_opt_bits &= ~SAVE_OPT_SLUUD;
                    }
                }
            }
            


#ifndef TARGET_API_LIB
            /*
            inchi_ios_eprint(stderr, "%ld: %s\r", num_inp, strHdr? strHdr : "");
            inchi_ios_eprint(pLog,  "%ld: %s\n", num_inp, strHdr? strHdr : "");
            */

            if ( !ip->bNoStructLabels && !(bInChI2Structure && (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)) ) 
            { 
                /* Added 2nd item: Do not output this extra line into the output SDfile. 2008-11-17 DCh */
                if ( strHdr && strstr(strHdr, "Structure:")) 
                {
                    inchi_ios_print(pOut,  "%s%s", strHdr, pTAB); /* output header */
#if ( FIX_DALKE_BUGS == 1 )
#else
                    sprintf(szMessage, "%s (%ld)",  strHdr? strHdr : "", num_inp);
#endif
                } 
                else 
                {
                    /*
                    inchi_ios_print(pOut,  "Structure %ld (%s)%s", num_inp, strHdr? strHdr : "", pTAB);
                    sprintf(szMessage, "Structure %ld (%s)%s", num_inp, strHdr? strHdr : "" , pTAB);
                    */
                    inchi_ios_print(pOut,  "Structure: %ld. (%s)%s", num_inp, strHdr? strHdr : "No struct name" , pTAB); /* output header */
#if ( FIX_DALKE_BUGS == 1 )
#else
                    sprintf(szMessage, "Structure: %ld. (%s)%s", num_inp, strHdr? strHdr : "No struct name", pTAB);
#endif
                }
                if ( strHdr && strHdr[0] ) 
                {
                    strncpy(ip->szSdfDataHeader, strHdr, sizeof(ip->szSdfDataHeader));
                    ip->szSdfDataHeader[sizeof(ip->szSdfDataHeader)-1] = '\0';
                    ip->pSdfLabel = NULL;
                    ip->pSdfValue = ip->szSdfDataHeader;
                } 
                else 
                {
                    ip->pSdfValue = NULL;
                    ip->szSdfDataHeader[0] = '\0';
                }
            }

#if ( FIX_DALKE_BUGS == 1 )
            sprintf(szMessage, "%ld: %.400s", num_inp, strHdr? strHdr : "");
#else
            sprintf(szMessage, "%ld: %s", num_inp, strHdr? strHdr : "");
#endif
#endif
            
            nInitLenMessage = strlen(szMessage);
            if ( strHdr ) 
            {
                szCurHdr = strHdr;
                strHdr = NULL;
            }
            if ( szCurHdr && ip && ip->first_struct_number > 0 ) 
            {
                /* check whether the structure should be skipped */
                static char szStruct[] = "Structure:";
                char *pStrNum = strstr(szCurHdr, szStruct);
                long cur_struct_number;
                if ( pStrNum ) 
                {
                    pStrNum += sizeof(szStruct)-1; /* -1 takes care of the string terminal zero */
                    cur_struct_number = inchi_strtol(pStrNum, NULL, 10);
                    if ( cur_struct_number ) 
                    {
                        OneInput.num_inp = cur_struct_number;
                    }
                    /* process request to bypass first several InChIs */
                    if ( cur_struct_number > 0 && cur_struct_number < ip->first_struct_number ) 
                    {

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )                   
                        inchi_fprintf(stderr, "Skipping %s\r", szMessage);
#endif
                        FreeInpInChI(&OneInput);
                        if ( szCurHdr ) 
                        {
                            inchi_free(szCurHdr);
                            szCurHdr = NULL;
                        }
                        INCHI_HEAPCHK
                        continue;
                    }
                }
            }
            
            num_processed ++;
            
            /* In case of splitting InChI into separate components */
            /* decide whether to keep /p in each component or      */
            /* output /p and /i/h as a separate component          */
            /* Note: if InChI is not to be splitted DO NOT create  */
            /* a separate component for /p, /i/h: it would be a bug*/
        
            InchiTimeGet(&ulTStart);

            INCHI_HEAPCHK
            ret = SetProtonsAndXchgIsoH(bInChI2Structure, bReqSplitOutputInChI, bReqProtonsForEachComponent,
                           bReqNonTaut, bReqStereo, num_components, nModeProtonIsoExchgH, &OneInput);
            INCHI_HEAPCHK
            if ( ret < 0 ) 
            {
                num_errors ++;
                goto exit_error;
            }

            sd->num_components[INCHI_BAS] = num_components[INCHI_BAS];
            sd->num_components[INCHI_REC] = num_components[INCHI_REC];

            /* do we have reconnected InChI ? */
            if ( (OneInput.nNumComponents[INCHI_REC][TAUT_YES] ||
                  OneInput.nNumComponents[INCHI_REC][TAUT_NON]) &&
                 (ip->bTautFlags & TG_FLAG_RECONNECT_COORD) &&
                 (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD) ) 
            {
                /* needed for InChI string output to include reconnected InChI */
                sd->bTautFlagsDone[0] |= TG_FLAG_DISCONNECT_COORD_DONE;
                bHasSomeReconnected = 1;
            }
            /* Do we have fixed H InChI ? */
            if ( bReqNonTaut &&
                 /*OneInput.nNumComponents[bHasSomeReconnected?INCHI_REC:INCHI_BAS][TAUT_NON]*/
                 (OneInput.nNumComponents[INCHI_REC][TAUT_NON] || 
                  OneInput.nNumComponents[INCHI_BAS][TAUT_NON]) ) 
            {
                bHasSomeFixedH = 1;
            }
            
            ulProcessingTime += InchiTimeElapsed(&ulTStart);


           
            if ( bInChI2Structure && !bInChI2InChI ) 
            {
                /**********************************************************************/
                /*                     InChi --> Structure                            */
                /**********************************************************************/
                

                int bINChIOutputOptions = 
#if ( I2S_MODIFY_OUTPUT == 1 )
                /* transfer user's InChI output options to serialization 10-12-2007 */
                    ip_inp->bINChIOutputOptions & (
                            INCHI_OUT_NO_AUX_INFO           |   /* do not output Aux Info */
                            INCHI_OUT_SHORT_AUX_INFO        |   /* output short version of Aux Info */
                            INCHI_OUT_ONLY_AUX_INFO         |   /* output only Aux Info */
                         /* INCHI_OUT_EMBED_REC             |*/   /* embed reconnected INChI into disconnected INChI */
                            INCHI_OUT_SDFILE_ONLY           |   /* save input data in a Molfile instead of creating INChI */
                            INCHI_OUT_XML                   |   /* output xml INChI */
                            INCHI_OUT_PLAIN_TEXT            |   /* output plain text INChI */
                            INCHI_OUT_PLAIN_TEXT_COMMENTS   |   /* output plain text annotation */
                            INCHI_OUT_XML_TEXT_COMMENTS     |   /* output xml text annotation */
                         /* INCHI_OUT_WINCHI_WINDOW         |*/   /* output into wINChI text window */
                            INCHI_OUT_TABBED_OUTPUT         |   /* tab-delimited (only for plain text) */
                            INCHI_OUT_SDFILE_ATOMS_DT       |   /* SDfile output H isotopes as D and T */
                            INCHI_OUT_SDFILE_SPLIT          |   /* Split SDfile into components */
                            0);

            
#else
                            0;
#endif


                SRM srm; /* rules how to handle bonds to metal atoms */
                StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM];


                /* prepare parameters */

                InchiTimeGet(&ulTStart);

                if ( bInputInStdFormat )
                {
                    if ( ip_inp->bINChIOutputOptions & INCHI_OUT_STDINCHI )
                        bINChIOutputOptions |= INCHI_OUT_STDINCHI;            
                }
                else
                {
                    if ( ip_inp->bINChIOutputOptions & INCHI_OUT_SAVEOPT ) 
                        bINChIOutputOptions |= INCHI_OUT_SAVEOPT;            
                }



                memset(pStruct, 0, sizeof(pStruct));
                
                /* structure restore parms */
                SetUpSrm(&srm);
                
                /* eliminate Fixed-H InChI that are exactly same as the corresponding Mobile-H structures */
                RemoveFixHInChIIdentical2MobH(&OneInput);

                /*-- recheck layers after the elimination; get optional stereo flags --*/
                ret = DetectInpInchiCreationOptions(&OneInput, 
                                                    &bHasSomeReconnected, &bHasMetal,
                                                    &bHasSomeFixedH, &nModeFlagsStereo, 
                                                    &bTautFlags);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, 
                              "Error in detecting input InChI options", "; ");
                    num_errors ++;
                    goto dealloc;
                }
                if ( bHasSomeFixedH && !bReqNonTaut ) 
                {
                    bHasSomeFixedH = 0;
                }
                /*---------------- set stereo flags ---------------------*/
                ip->nMode &= ~(REQ_MODE_STEREO | REQ_MODE_ISO_STEREO |
                               REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO |
                               REQ_MODE_CHIR_FLG_STEREO |
                               REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
                ip->nMode |= nModeFlagsStereo;
                
                /* Remove Phosphine and Arsine Stereo Flags */
                ip->bTautFlags &= ~TG_FLAG_PHOSPHINE_STEREO;
                ip->bTautFlags &= ~TG_FLAG_ARSINE_STEREO;
                ip->bTautFlags &= ~TG_FLAG_FIX_SP3_BUG;

                ip->bTautFlags |= bTautFlags;

                /* mark Disconnected InChI components that are exactly came as Reconnected ones */
                /* Disconnected will have a negative number of the reconnected component */
                /* Reconnected will have a positive number of the disconnected component */
                MarkDisconectedIdenticalToReconnected (&OneInput);

                /*****************************************************************************/
                /* pay attention to:                                                         */
                /* 1) .nLink < 0 in Disonnected which means InChI is same as in Reconnected  */
                /*    The component in Reconnected has .nLink pointing to the Disconnected;  */
                /*    each .nLink = (1+component index) or -(1+component index)              */
                /*    In the future .nLink>0 in Disconnected shall point to the Reconnectrd  */
                /*    component from which it was created                                    */
                /* 2) Currently reversed structures from Disconnected components are created */
                /*    and abandoned if Reconnected layer exists                              */
                /* 3) Connect/disconnect H depends on the presence of atom/bond parity       */
                /*    The combined Mobile/Fixed-H parity should be set for Fixed-H components*/
                /* 4) No comparison of the Disconnected layer is done if Reconnected exists  */
                /* 5) Reading InChI was not fully tested in case one component has stereo in */
                /*    both Mobile-H and Fixed-H layers while another component has stereo    */
                /*    only in Mobile-H layer                                                 */
                /*****************************************************************************/

                /* main conversion InChI->Structure for each component and  */
                /* after that pStruct[iRec][iMobH][iComponent].at2 is the structure,      */
                /* pStruct[iRec][iMobH][iComponent].RevInChI full InChI for the structure */
                /* In case of both Fixed-H and Mobile-H layers the results are in iMobH=0 */
                /* In case of only Mobile-H/Main layer the results are in iMobH=1         */
                ulProcessingTime += InchiTimeElapsed(&ulTStart);

                sd->ulStructTime = 0;
                ret = AllInchiToStructure(ip, sd, num_inp, szCurHdr, &srm, bHasSomeFixedH, pStruct, &OneInput);

                ulProcessingTime += sd->ulStructTime;
                InchiTimeGet(&ulTStart);
                /* ret < 0 is error code; ret > 0 is number of errors */
                /* in pStruct[iInchiRec][iMobileH][iComponent].nError */

                if ( ret) 
                {
                    /* conversion error */
                    num_errors ++;
                    goto dealloc;
                }
                

                /* an attempt to fix the numumber of removed protons in case of Mobile-H */
                if ( !OneInput.nNumProtons[INCHI_BAS][TAUT_YES].pNumProtons &&
                     !OneInput.nNumProtons[INCHI_REC][TAUT_YES].pNumProtons ) 
                {
                    ret = AddProtonAndIsoHBalanceToMobHStruct(ip, sd, num_inp, bHasSomeFixedH, szCurHdr, pStruct, &OneInput);
                    if ( ret < 0 ) 
                    {
                        AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Add/Remove protons error", "; ");
                        num_errors ++;
                        goto dealloc;
                    }
                }
                
                /* compare InChI from the Reversed Structure to the original input InChI */
                
                ret = CompareAllOrigInchiToRevInChI(pStruct, &OneInput, bHasSomeFixedH, num_inp, szCurHdr);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "InChI compare error", "; ");
                    num_errors ++;
                    goto dealloc;
                }
                
                ret = CompareAllDisconnectedOrigInchiToRevInChI(pStruct, &OneInput, bHasSomeFixedH,
                                                                 num_inp, szCurHdr);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "InChI compare2 error", "; ");
                    num_errors ++;
                    goto dealloc;
                }

                if ( WarningFlags ) 
                {
                    for ( i = 0; i < 2; i ++ ) 
                    {
                        for ( j = 0; j < TAUT_NUM; j ++ ) 
                        {
                            WarningFlags[i][j] = (unsigned long)OneInput.CompareInchiFlags[i][j];
                        }
                    }
                }
                
                ulProcessingTime += InchiTimeElapsed(&ulTStart);

#ifndef COMPILE_ANSI_ONLY
                ret = DisplayStructureComponents(ip, sd, num_inp, szCurHdr,
                                                  &srm, bReqNonTaut, pStruct, &OneInput);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Display structure error", "; ");
                }
#endif
                InchiTimeGet(&ulTStart);
                ret = MergeStructureComponents(ip, sd, num_inp, szCurHdr, &srm, bReqNonTaut, pStruct, &OneInput);
                ulProcessingTime += InchiTimeElapsed(&ulTStart);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Merge Components error", "; ");
                    num_errors ++;
                    goto dealloc;
                }



#ifdef TARGET_API_LIB
/*------------- for debug only -------------------
                InchiTimeGet(&ulTStart);
                ret = OutputInChIOutOfStrFromINChI(ip, sd, num_inp, 0, 
                                                   pOut, pLog, &OneInput,
                                                   save_opt_bits);
                ulProcessingTime += InchiTimeElapsed(&ulTStart);
                if ( ret < 0 ) {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Restored structure to InChI conversion failed", "; ");
                    goto dealloc;
                }
-------------------------------------------------*/
                if ( at && num_at ) 
                {
                    *at     = OneInput.atom;
                    *num_at = OneInput.num_atoms;
                    OneInput.atom = NULL;
                }
#else
                InchiTimeGet(&ulTStart);
                ret = OutputInChIOutOfStrFromINChI(ip, sd, num_inp, bINChIOutputOptions, 
                                                   pOut, /*pLog*/ NULL, &OneInput, 
                                                   bHasSomeFixedH, save_opt_bits);
                ulProcessingTime += InchiTimeElapsed(&ulTStart);
                if ( ret < 0 ) 
                {
                    AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Restored structure to InChI conversion error", "; ");
                    num_errors ++;
                    goto dealloc;
                }
#endif

                if ( 1 ) 
                {
                    int len;
                    InchiTimeGet(&ulTStart);
                    FillOutCompareMessage(szMessage, nMessageLen, OneInput.CompareInchiFlags[0]);
                    if ( OneInput.CompareInchiFlags[1][0] || OneInput.CompareInchiFlags[1][1] ) 
                    {
                        AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "Disconnected: ", "; ");
                        FillOutCompareMessage(szMessage, nMessageLen, OneInput.CompareInchiFlags[1]);
                    }
                    /* add a metal warning */
                    if ( bHasMetal && nInitLenMessage < (len=(int)strlen(szMessage)) ) 
                    {
                        char szMetal[] = " (Metal compound)";
                        int shift;
                        if ( len + (int)sizeof(szMetal) > nMessageLen ) {
                            len = nMessageLen - (int)sizeof(szMetal);
                        }
                        shift = nInitLenMessage + (int)sizeof(szMetal) - 1;
                        memmove(szMessage+shift, szMessage + nInitLenMessage, (len-nInitLenMessage)*sizeof(szMessage[0]));
                        memcpy(szMessage + nInitLenMessage, szMetal, sizeof(szMetal)-sizeof(szMessage[0]));
                        szMessage[shift+len-nInitLenMessage] = '\0';
                    }
                    ulProcessingTime += InchiTimeElapsed(&ulTStart);
                }
                
                ret = 0;

                /* deallocate */
dealloc:        
                if ( ret ) 
                {
                    if ( ret > 0 ) 
                    {
                        int iRec, iMob, iComp, nComp, len;
                        char szTemp[128];
                        AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "*Conversion failed on component(s)", "; ");
                        len = strlen(szMessage);
                        for ( iRec = 0; iRec < INCHI_NUM; iRec ++ ) 
                        {
                            for ( iMob = bHasSomeFixedH? TAUT_NON : TAUT_YES; iMob < TAUT_NUM; iMob ++ ) 
                            {
                                nComp = OneInput.nNumComponents[iRec][iMob];
                                if ( !pStruct[iRec][iMob] ) 
                                {
                                    continue;
                                }
                                for ( iComp = 0; iComp < nComp; iComp ++ ) 
                                {
                                    if ( pStruct[iRec][iMob][iComp].nError ) 
                                    {
                                        char *szFormula = OneInput.pInpInChI[iRec][iMob][iComp].szHillFormula;
                                        sprintf (szTemp,
#if ( FIX_DALKE_BUGS == 1 )                                            
                                            " %s%s%d(%.96s)",
#else
                                            " %s%s%d(%s)",
#endif
                                            !bHasSomeReconnected? "" : iRec? "R" : "D", 
                                            !bHasSomeFixedH? "": iMob? "M" : "F",
                                            iComp + 1, szFormula? szFormula : "???");
                                        AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, szTemp, NULL);
                                    }
                                }
                            }
                        }
                    } 
                    else 
                    {
                        if ( ret == CT_USER_QUIT_ERR ) 
                        {
                            AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "*Terminated by the user*", "; ");
                        } else {
                            AddOneMsg(szMessage, (int)strlen(szMessage), nMessageLen, "*Conversion failed*", "; ");
                        }
                    }
                }
                InchiTimeGet(&ulTStart);
                /* print one structure report */
                if ( szMsg && nMsgLen > 1 ) 
                {
                    int len = inchi_min( (int)strlen(szMessage), nMsgLen-1);
                    if ( len > 0 ) 
                    {
                        memcpy( szMsg, szMessage, len);
                        szMsg[len] = '\0';
                    } 
                    else 
                    {
                        szMsg[0] = '\0';
                    }
                }
                if ( nInitLenMessage < (int)strlen(szMessage) ) 
                {
                    inchi_ios_eprint(pLog, "%s\n", szMessage);
                }
#ifndef TARGET_API_LIB
                else 
                {
                    /*^^^inchi_ios_eprint( stderr, "%s\r", szMessage );*/
                    inchi_fprintf( stderr, "%s\r", szMessage );
                }
#endif                


                FreeStrFromINChI( pStruct, OneInput.nNumComponents );
                FreeInpInChI( &OneInput );
                if ( szCurHdr ) 
                {
                    inchi_free( szCurHdr );
                    szCurHdr = NULL;
                }
                INCHI_HEAPCHK
                ulProcessingTime += InchiTimeElapsed( &ulTStart );
                if ( ret < 0 ) 
                {
                    goto exit_error;
                }

            }  /* if ( bInChI2Structure && !bInChI2InChI ) */
            


            else if ( !bInChI2Structure && bInChI2InChI ) 
            {            
                /**********************************************************************/
                /*                     InChi --> InChI string(s)                      */
                /**********************************************************************/

                int tmp = ip->bNoStructLabels;
                
                InchiTimeGet( &ulTStart );
                
                ip->bNoStructLabels = 1;
                INCHI_HEAPCHK
                ip->pSdfValue = NULL;
                ip->pSdfLabel = NULL;
#if ( FIX_DALKE_BUGS == 1 )
                SetHillFormFromInChI( &OneInput );
#endif


                ret = OutputInChIAsRequested(pOut, pLog, ip, sd, &OneInput, 
                                             num_components, nModeProtonIsoExchgH, 
                                             num_inp, save_opt_bits);


#if ( !defined(TARGET_API_LIB) && defined(TARGET_EXE_STANDALONE) )
                /*^^^ calculate InChIKey if requested */
                /* However, do not calculat/write it if this function is called from within dll */
                {
                    char ik_string[256];    /*^^^ Resulting InChIKey string */
                    int ik_ret=0;           /*^^^ InChIKey-calc result code */
                    int xhash1, xhash2;
                    char szXtra1[65], szXtra2[65];

                    inchi_ios_flush2(pLog, stderr);

                    /*^^^ post-1.02b addition - correctly treat tabbed output with InChIKey */
                    if ( ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT ) 
                        if ( ip->bCalcInChIHash != INCHIHASH_NONE )
                            if (pOut->s.pStr)
                                if (pOut->s.nUsedLength>0)
                                    if (pOut->s.pStr[pOut->s.nUsedLength-1]=='\n')
                                        /* replace LF with TAB */
                                        pOut->s.pStr[pOut->s.nUsedLength-1] = '\t';


                    if ( ip->bCalcInChIHash != INCHIHASH_NONE )
                    {
                        char *buf = NULL;
                        size_t slen = pOut->s.nUsedLength;
                        extract_inchi_substring(&buf, pOut->s.pStr, slen);
            
                        if (NULL!=buf)
                        {
                            xhash1 = xhash2 = 0;
                            if ( ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1 ) ||
                                 ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ) )
                                 xhash1 = 1;
                            if ( ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2 ) ||
                                 ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ) )
                                 xhash2 = 1;                
                            ik_ret = GetINCHIKeyFromINCHI(buf, xhash1, xhash2, 
                                                          ik_string, szXtra1, szXtra2);
                            inchi_free(buf);
                        }
                        else
                            ik_ret = INCHIKEY_NOT_ENOUGH_MEMORY;                     
        

                        if (ik_ret==INCHIKEY_OK)   
                        {
                            inchi_ios_print(pOut, "InChIKey=%-s\n",ik_string);
                        }
                        else    
                        {
                            inchi_ios_print(pLog, "Warning (Could not compute InChIKey: ", num_inp);
                            switch(ik_ret)
                            {
                                case INCHIKEY_UNKNOWN_ERROR:
                                    inchi_ios_print(pLog, "unresolved error)");
                                    break;
                                case INCHIKEY_EMPTY_INPUT:
                                    inchi_ios_print(pLog,  "got an empty string)");
                                    break;
                                case INCHIKEY_INVALID_INCHI_PREFIX:
                                case INCHIKEY_INVALID_INCHI:
                                case INCHIKEY_INVALID_STD_INCHI:
                                    inchi_ios_print(pLog, "got non-InChI string)");
                                    break;
                                case INCHIKEY_NOT_ENOUGH_MEMORY:
                                    inchi_ios_print(pLog, "not enough memory to treat the string)");
                                    break;
                                default:inchi_ios_print(pLog, "internal program error)");
                                    break;
                            }
                            inchi_ios_print(pLog, " structure #%-lu.\n", num_inp);
                            if ( ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT )
                                inchi_ios_print(pOut, "\n");
                        } /* if (ip->bCalcInChIHash!=INCHIHASH_NONE) */

                
                        inchi_ios_flush(pOut);
                        inchi_ios_flush2(pLog, stderr);
            
                    } 

                    else
                        inchi_ios_flush(pOut);
                } /* calculate InChIKey if requested */
#endif

                ip->bNoStructLabels = tmp;

#ifndef TARGET_API_LIB
                if ( ret < 0 ) 
                {
                    if ( szCurHdr && szCurHdr[0] ) 
                    {
                        inchi_ios_eprint( pLog, "Error %d creating InChI string %s\n", ret, szCurHdr );
                    } 
                    else 
                    {
                        inchi_ios_eprint( pLog, "Error %d creating InChI string, Structure %ld\n", ret, num_inp );
                    }
                    num_errors ++;
                } 
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )                   
                else
                if ( szCurHdr && szCurHdr[0] ) 
                {
                    inchi_fprintf( stderr, "%s\r", szCurHdr );
                }
#endif
#endif

                if ( szCurHdr ) 
                {
                    inchi_free( szCurHdr );
                    szCurHdr = NULL;
                }
        
                INCHI_HEAPCHK
                ulProcessingTime += InchiTimeElapsed( &ulTStart );
        
            } /* if ( !bInChI2Structure && bInChI2InChI )  */
        
            else 
            {
                inchi_ios_eprint( pLog, "\nWrong command line options: expected Inch2Struct or Inchi2Inchi\n", num_inp );
                break;
            }
    

            if ( nReadStatus == RI_ERR_EOF ) 
            {
                break;
            }

        
        } /* InChI has been successfully read */
    
        else 
    
        {
    
            /* InChI could not be read */
            if ( nReadStatus == RI_ERR_EOF && nErr == 0 && pState == 0 && !strHdr ) 
            {
                inchi_ios_eprint( pLog, "\nEnd of file detected after structure %ld.    \n", num_inp );
            } 
            else 
            {
                /* output InChI parsing error message */
                char szHdrSimulation[128];
                num_inp ++;
                sprintf( szHdrSimulation, "Structure: %ld", num_inp );
                inchi_ios_eprint( pLog, "\n%s %s (%d) in %s (%d)\n", strHdr? strHdr : szHdrSimulation, getInchiErrName(nErr), nErr,  getInchiStateReadErr(pState), pState );
                num_errors ++;
                num_processed ++;
            }
            if ( strHdr ) 
            {
                inchi_free( strHdr );
                strHdr = NULL;
            }
            if ( szCurHdr ) 
            {
                inchi_free( szCurHdr );
                szCurHdr = NULL;
            }
            FreeInpInChI( &OneInput );
        }


#ifdef TARGET_EXE_STANDALONE
#ifndef TARGET_API_LIB
        inchi_ios_flush(pOut);
        inchi_ios_flush2(pLog, stderr);
#endif
#endif


#ifdef TARGET_API_LIB
        break;  /* exit after the 1st structure */
#endif

    
    } /* while */




exit_error:

    FreeInpInChI( &OneInput );
    if ( strHdr ) 
    {
        inchi_free( strHdr );
        strHdr = NULL;
    }
    if ( pLine->str ) 
    {
        inchi_free( pLine->str );
    }
    if ( szCurHdr ) 
    {
        inchi_free( szCurHdr );
        szCurHdr = NULL;
    }
    
    INCHI_HEAPCHK

    if ( sd_inp ) 
    {
        sd_inp->ulStructTime = ulProcessingTime;
        sd_inp->fPtrStart    = num_processed;
        sd_inp->fPtrEnd      = num_errors;
    }
    return ret;
}



/**********************************************************************************************/
int OutputInChIAsRequested(INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, 
                           ICHICONST INPUT_PARMS *ip_inp, 
                           STRUCT_DATA *sd_inp, 
                           InpInChI *OneInput, 
                           int num_components[INCHI_NUM],
                           MODE_PIXH nModeProtonIsoExchgH[INCHI_NUM], 
                           long num_inp, unsigned char save_opt_bits)
{
    int      j, k, k1, k2, ret2=0, iINChI, iINChI1, iINChI2;
    PINChI2 *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];
    int bReqNonTaut;
    int bHasSomeReconnected;

    INPUT_PARMS ip_local;
    STRUCT_DATA sd_local;
    INPUT_PARMS *ip = &ip_local;
    STRUCT_DATA *sd = &sd_local;
    NORM_CANON_FLAGS ncFlags;
    NORM_CANON_FLAGS *pncFlags = &ncFlags;
    const int nStrLen = NSTRLEN;
    char *pStr = NULL;
    int  nRet1, bSortPrintINChIFlags;
    int  bReqSplitOutputInChI;
    int  nNumOutputComponents;
    
    nRet1   = 0;
    k1 = k2 = 0;
    memset( pncFlags, 0, sizeof(*pncFlags) );
    memset( pINChI, 0, sizeof(pINChI) );
    memset( pINChI_Aux, 0, sizeof(pINChI_Aux) );

    *ip = *ip_inp;
    *sd = *sd_inp;
    bHasSomeReconnected = 0;
    bSortPrintINChIFlags = 0;
    nNumOutputComponents = 0;
    bReqNonTaut = (0 != (ip->nMode & REQ_MODE_BASIC));
    bReqSplitOutputInChI = (0 != (ip->bReadInChIOptions & READ_INCHI_SPLIT_OUTPUT));

    INCHI_HEAPCHK

    if ( num_components[INCHI_BAS] ) 
    {
        MYREALLOC2(PINChI2, PINChI_Aux2, pINChI[INCHI_BAS], pINChI_Aux[INCHI_BAS], num_components[INCHI_BAS], num_components[INCHI_BAS], k1);
    }
    if ( num_components[INCHI_REC] ) 
    {
        MYREALLOC2(PINChI2, PINChI_Aux2, pINChI[INCHI_REC], pINChI_Aux[INCHI_REC], num_components[INCHI_REC], num_components[INCHI_REC], k2);
    }
    pStr = (char*) inchi_malloc( nStrLen * sizeof(pStr[0]) );

    INCHI_HEAPCHK

    if ( k1 || k2 || !pStr ) 
    {
        ret2 = RI_ERR_ALLOC;
        goto exit_error;
    }
    
    if ( num_components[INCHI_REC] &&
         (ip->bTautFlags & TG_FLAG_RECONNECT_COORD) &&
         (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD) ) 
    {
        sd->bTautFlagsDone[0] |= TG_FLAG_DISCONNECT_COORD_DONE;
        bHasSomeReconnected = 1;
    }


    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) 
    {
        for ( j = 0; j < TAUT_NUM; j ++ ) 
        {
            if ( bReqNonTaut || (j != TAUT_NON && OneInput->pInpInChI[iINChI][j]) )
            {
                for ( k = 0; k < num_components[iINChI]; k ++ ) 
                {
                    /* allocate InChI & AuxInfo */
                    if ( !(pINChI[iINChI][k][j] = (INChI *) inchi_calloc(1, sizeof(INChI)) ) ) 
                    {
                        ret2 = RI_ERR_ALLOC;
                        goto exit_error;
                    }
                    if ( !(pINChI_Aux[iINChI][k][j] = (INChI_Aux *) inchi_calloc(1, sizeof(INChI_Aux)) ) ) 
                    {
                        ret2 = RI_ERR_ALLOC;
                        goto exit_error;
                    }
                    /* copy InChI & AuxInfo */
                    if ( k < OneInput->nNumComponents[iINChI][j] ) 
                    {

                        /* copy InChI */
                        *pINChI[iINChI][k][j] = OneInput->pInpInChI[iINChI][j][k];
                        memset(&OneInput->pInpInChI[iINChI][j][k], 0, sizeof(OneInput->pInpInChI[iINChI][j][k]));
                        INCHI_HEAPCHK
                        /* take care of protons in AuxInfo */

                        if ( nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_EACH && j == TAUT_YES ) 
                        {
                            pINChI_Aux[iINChI][k][j]->nNumRemovedProtons = 
                                   OneInput->nNumProtons[iINChI][j].pNumProtons[k].nNumRemovedProtons;
                            for ( k1 = 0; k1 < NUM_H_ISOTOPES; k1 ++ ) 
                            {
                                pINChI_Aux[iINChI][k][j]->nNumRemovedIsotopicH[k1] =
                                    OneInput->nNumProtons[iINChI][j].pNumProtons[k].nNumRemovedIsotopicH[k1];
                            }
                            INCHI_HEAPCHK
                        } 
                        else if ( (!k && nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_FIRST) ||
                             (k+1 == OneInput->nNumComponents[iINChI][j] &&
                             nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_A_PIXH_COMPONENT) ) 
                        {
                            /* add protons and exchangeable isotopic H to the first component's AuxInfo */
                            pINChI_Aux[iINChI][k][j]->nNumRemovedProtons = OneInput->nNumProtons[iINChI][j].nNumRemovedProtons;
                            for ( k1 = 0; k1 < NUM_H_ISOTOPES; k1 ++ ) 
                            {
                                pINChI_Aux[iINChI][k][j]->nNumRemovedIsotopicH[k1] =
                                    OneInput->nNumProtons[iINChI][j].nNumRemovedIsotopicH[k1];
                            }
                            INCHI_HEAPCHK
                        } 
                        else 
                        {
                            pINChI_Aux[iINChI][k][j]->bDeleted = pINChI[iINChI][k][j]->bDeleted;
                        }

                        if ( j == TAUT_YES && pINChI[iINChI][k][j] && pINChI[iINChI][k][j]->nNumberOfAtoms &&
                             !pINChI[iINChI][k][j]->nNum_H_fixed ) 
                        {
                            /* serializer crashes if it is not allocated */
                            pINChI[iINChI][k][j]->nNum_H_fixed = (S_CHAR *)inchi_calloc(pINChI[iINChI][k][j]->nNumberOfAtoms+1, sizeof(pINChI[0][0][0]->nNum_H_fixed[0]) );
                        }

                        if ( j == TAUT_YES && k < OneInput->nNumComponents[iINChI][TAUT_NON] &&
                             pINChI[iINChI][k][j] && pINChI[iINChI][k][j]->nNumberOfAtoms &&
                             pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms &&
                             !CompareReversedINChI( pINChI[iINChI][k][j], pINChI[iINChI][k][TAUT_NON], NULL, NULL ) ) {
                            pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms = 0; /* eliminate non-taut equal to taut */
                        }
                             

                    } 
                    else 
                    {
                        /* extra component, usually it is a Mobile H component */
                        /* corresponding to a free proton component in Fixed H */
                        pINChI[iINChI][k][j]->bDeleted = 1;
                        pINChI_Aux[iINChI][k][j]->bDeleted = 1;
                    }

                } /* k */

            } /* if ( bReqNonTaut || j != TAUT_NON && OneInput->pInpInChI[iINChI][j] )  */

            if ( OneInput->pInpInChI[iINChI][j] ) 
            {
                INCHI_HEAPCHK
                inchi_free(OneInput->pInpInChI[iINChI][j]);
                OneInput->pInpInChI[iINChI][j] = NULL;
            }
        
        } /* j */

    } /* iINChI */
    
    if ( bReqSplitOutputInChI ) 
    {
        if ( bHasSomeReconnected ) 
        {
            iINChI1 = INCHI_REC; /* only reconnected */
            iINChI2 = INCHI_NUM;
            sd->num_components[INCHI_BAS] = sd->num_components[INCHI_REC];
        } 
        else 
        {
            iINChI1 = 0;         /* only disconnected */
            iINChI2 = iINChI1+1;
        }
        sd->num_components[INCHI_REC] = 0;  /* treat reconnected as connected */
        nNumOutputComponents = sd->num_components[INCHI_BAS];
    } 
    else 
    {
        iINChI1 = 0;
        iINChI2 = INCHI_NUM;
        nNumOutputComponents = 1;
    }


    for ( k1 = 0, k2 = (bReqSplitOutputInChI? k1+1 : nNumOutputComponents); k1 < k2 && k1 < nNumOutputComponents; k1=k2, k2 ++ ) 
    {

        if ( bReqSplitOutputInChI ) 
        {
            sd->num_components[INCHI_BAS] = 1;
            sd->num_components[INCHI_REC] = 0;
            /* additional data */
            sd->num_non_taut[INCHI_BAS] =
            sd->num_taut[INCHI_BAS]     =
            sd->num_non_taut[INCHI_REC] =
            sd->num_taut[INCHI_REC]     = 0;
            iINChI = iINChI1;
            for ( j = 0; j < TAUT_NUM && sd->num_components[iINChI]; j ++ ) 
            {
                for ( k = k1; k < k2; k ++ ) 
                {
                    /*  find where the current processed structure is located */
                    int cur_is_in_non_taut = (pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms>0);
                    int cur_is_in_taut     = (pINChI[iINChI][k][TAUT_YES] && pINChI[iINChI][k][TAUT_YES]->nNumberOfAtoms>0);
                    int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[iINChI][k][TAUT_NON]->lenTautomer) ||
                                          (cur_is_in_taut     && 0 == pINChI[iINChI][k][TAUT_YES]->lenTautomer);
                    int cur_is_taut     = cur_is_in_taut     && 0 <  pINChI[iINChI][k][TAUT_YES]->lenTautomer;
                    if ( cur_is_non_taut + cur_is_taut ) 
                    {
                        /*  count tautomeric and non-tautomeric components of the structures */
                        /*
                        int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
                        int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
                        */
                        sd->num_non_taut[INCHI_BAS] += cur_is_non_taut;
                        sd->num_taut[INCHI_BAS]     += cur_is_taut;
                    }
                }
            }
            INCHI_HEAPCHK
        } 
        else 
        {
            sd->num_components[INCHI_BAS] = inchi_max( OneInput->nNumComponents[INCHI_BAS][TAUT_YES],
                                                       OneInput->nNumComponents[INCHI_BAS][TAUT_NON] );
            sd->num_components[INCHI_REC] = inchi_max( OneInput->nNumComponents[INCHI_REC][TAUT_YES],
                                                       OneInput->nNumComponents[INCHI_REC][TAUT_NON] );
        
            /* additional data needed for SortAndPrintINChI() */
            for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) 
            {
                sd->num_non_taut[iINChI] =
                sd->num_taut[iINChI]     = 0;
                for ( j = 0; j < TAUT_NUM && sd->num_components[iINChI]; j ++ ) 
                {
                    for ( k = k1; k < k2; k ++ ) 
                    {
                        /*  find where the current processed structure is located */
                        int cur_is_in_non_taut = (pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms>0);
                        int cur_is_in_taut     = (pINChI[iINChI][k][TAUT_YES] && pINChI[iINChI][k][TAUT_YES]->nNumberOfAtoms>0);
                        int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[iINChI][k][TAUT_NON]->lenTautomer) ||
                                              (cur_is_in_taut     && 0 == pINChI[iINChI][k][TAUT_YES]->lenTautomer);
                        int cur_is_taut     = cur_is_in_taut     && 0 <  pINChI[iINChI][k][TAUT_YES]->lenTautomer;
                        if ( cur_is_non_taut + cur_is_taut ) {
                            /*  count tautomeric and non-tautomeric components of the structures */
                            /*
                            int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
                            int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
                            */
                            sd->num_non_taut[iINChI] += cur_is_non_taut;
                            sd->num_taut[iINChI]     += cur_is_taut;
                        }
                    }
                }
            }
            INCHI_HEAPCHK
        }
        if ( bReqSplitOutputInChI ) 
        {
            /* output components one by one (for splitting input InChI into components) */
            PINChI2 *pInChI_2[INCHI_NUM];
            PINChI_Aux2 *pInChI_Aux_2[INCHI_NUM];
            INChI *pInChI_1[1][2];
            INChI_Aux *pInChI_Aux_1[1][2];
            memset( pInChI_2, 0, sizeof(pInChI_2) );
            memset( pInChI_Aux_2, 0, sizeof(pInChI_Aux_2) );

            for ( j = 0; j < TAUT_NUM; j ++ ) 
            {
                pInChI_1[0][j]     = pINChI[iINChI1][k1][j];
                pInChI_Aux_1[0][j] = pINChI_Aux[iINChI1][k1][j];
            }
            pInChI_2[INCHI_BAS]     = pInChI_1;
            pInChI_Aux_2[INCHI_BAS] = pInChI_Aux_1;
            /* make sure purely reconnected InChI is marked as ReChI, not InChI */
            if ( bHasSomeReconnected &&
                 (bInChIHasReconnectedMetal( pInChI_1[0][TAUT_YES] ) ||
                  bInChIHasReconnectedMetal( pInChI_1[0][TAUT_NON] ) ) ) 
            {
                bSortPrintINChIFlags = FLAG_SORT_PRINT_ReChI_PREFIX;
            } 
            else 
            {
                bSortPrintINChIFlags = 0;
            }
            INCHI_HEAPCHK
            nRet1 = SortAndPrintINChI(pOut, pStr, nStrLen, pLog, ip, NULL /*orig_inp_data*/, NULL  /*prep_inp_data*/,
                                      NULL /*composite_norm_data*/, NULL /*pOrigStruct*/,
                                      sd->num_components, sd->num_non_taut, sd->num_taut,
                                      sd->bTautFlags, sd->bTautFlagsDone, pncFlags, num_inp,
                                      pInChI_2, pInChI_Aux_2, &bSortPrintINChIFlags,
                                      save_opt_bits);
            INCHI_HEAPCHK
        } 
        else 
        {
            INCHI_HEAPCHK
            bSortPrintINChIFlags = 0;
            nRet1 = SortAndPrintINChI(pOut, pStr, nStrLen, pLog, ip, NULL /*orig_inp_data*/, NULL  /*prep_inp_data*/,
                                      NULL /*composite_norm_data*/, NULL /*pOrigStruct*/,
                                      sd->num_components, sd->num_non_taut, sd->num_taut,
                                      sd->bTautFlags, sd->bTautFlagsDone, pncFlags, num_inp,
                                      pINChI, pINChI_Aux, &bSortPrintINChIFlags,
                                      save_opt_bits);
            INCHI_HEAPCHK
        }
        if ( nRet1 == _IS_FATAL || nRet1 == _IS_ERROR ) 
        {
            break;
        }
    }


    INCHI_HEAPCHK
    FreeAllINChIArrays( pINChI, pINChI_Aux, num_components );
    INCHI_HEAPCHK
    
    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) 
    {
        for ( j = 0; j < TAUT_NUM; j ++ ) 
        {
            if ( OneInput->nNumProtons[iINChI][j].pNumProtons ) 
            {
                inchi_free( OneInput->nNumProtons[iINChI][j].pNumProtons );
                OneInput->nNumProtons[iINChI][j].pNumProtons = NULL;
            }
        }
    }

    INCHI_HEAPCHK
    if ( nRet1 == _IS_FATAL || nRet1 == _IS_ERROR ) 
    {
        ret2 = RI_ERR_PROGR;
    }

exit_error:
    if ( pStr ) 
        inchi_free( pStr );
    return ret2;
}

/**************************************************************************************/
int GetNumNeighborsFromInchi( INChI *pInChI, AT_NUMB nAtNumber )
{
    int i, j, n_vertex, n_neigh, nNumNeigh, bTautAtom, nNumH, nTotNumNeigh, num_atoms;
    AT_NUMB  taut_at_number;
    nAtNumber -= 1;
    nNumNeigh = 0; /* number of bonds */
    bTautAtom = 0; /* 1 if atom belongs to a Mobile-H group */
    nNumH     = 0; /* number of terminal neighbors H */
    num_atoms = pInChI->nNumberOfAtoms;
    /* from RestoreAtomConnectionsSetStereo() */
    /* Connection table structure: 
       Vert(1) [, Neigh(11), Neigh(12),...], Vert(2) [, Neigh(2,1), Neigh(2,2),...] ...
       where Neigh(i,1) < Neigh(i,2) <... < Vert(i);
             Vert(i) < Vert(i+1)
     */
    for ( i = 1, n_vertex = pInChI->nConnTable[0]-1; i < pInChI->lenConnTable; i ++ ) {
        if ( (n_neigh = pInChI->nConnTable[i]-1) < n_vertex ) {
            /*  vertex - neighbor connection */
            nNumNeigh += ( nAtNumber == n_vertex || nAtNumber == n_neigh );
        } else
        /* n_neigh is the next vertex */
        if ( (n_vertex = n_neigh) >= num_atoms ) {
            return  RI_ERR_PROGR;
        }
    }
    /* is atom tautomeric, from GetTgroupInfoFromInChI() */
    if ( pInChI && pInChI->lenTautomer > 1 && pInChI->nTautomer && pInChI->nTautomer[0] > 0 ) {
        int itg, len_tg;
        int tot_len_tg = pInChI->lenTautomer - T_GROUP_HDR_LEN*pInChI->nTautomer[0] - 1; /* number of endpoints */
        j = 1; /* index in pInChI->nTautomer[] */
        i = 0; /* index in ti->nEndpointAtomNumber[] */
        for ( itg = 0; itg < pInChI->nTautomer[0]; itg ++ ) {
            len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
            j      += T_GROUP_HDR_LEN;   /* skip t-group header */
            len_tg -= T_GROUP_HDR_LEN-1;
            for( ; 0 < len_tg --; j ++, i ++ ) {
                taut_at_number = pInChI->nTautomer[j]-1; /* Mobile-H group atom number */
                bTautAtom += (taut_at_number == nAtNumber);
            }
        }
        if ( i != tot_len_tg ) {
            return RI_ERR_PROGR;
        }
    }
    /* count hydrogen neighbors */
    if ( pInChI->nNum_H ) {
        nNumH = pInChI->nNum_H[nAtNumber];
    }
    /* conclusion: if not tautomeric then return positive number, otherwise add 1000 */
    nTotNumNeigh = nNumNeigh + nNumH;
    if ( bTautAtom ) {
        nTotNumNeigh += 1000;
    }
    return nTotNumNeigh;
}
/**************************************************************************************/
int CountStereoTypes( INChI *pInChI, int *num_known_SB, int *num_known_SC,
                                     int *num_unk_und_SB, int *num_unk_und_SC,
                                     int *num_SC_PIII, int *num_SC_AsIII)
{
    static U_CHAR el_number_P=0, el_number_As=0;
    INChI_Stereo *Stereo;
    int           i, ret;
    AT_NUMB       nAtNumber;
    U_CHAR        el_number;
    if ( !pInChI->nNumberOfAtoms || pInChI->bDeleted ) {
        return 0; /* no InChI */
    }
    Stereo = (pInChI->StereoIsotopic &&
         (pInChI->StereoIsotopic->nNumberOfStereoBonds +
         pInChI->StereoIsotopic->nNumberOfStereoCenters ))? pInChI->StereoIsotopic:
                           (pInChI->Stereo &&
         (pInChI->Stereo->nNumberOfStereoBonds +
         pInChI->Stereo->nNumberOfStereoCenters ))? pInChI->Stereo : NULL;
    if ( !Stereo ) {
        return 1; /* No Stereo */
    }
    /* one-time initialization */
    if ( !el_number_P ) {
        el_number_P  = (U_CHAR)get_periodic_table_number( "P" );
        el_number_As = (U_CHAR)get_periodic_table_number( "As" );
    }
    /* count SB and cumulenes */
    for ( i = 0; i < Stereo->nNumberOfStereoBonds; i ++ ) {
        if ( ATOM_PARITY_WELL_DEF(Stereo->b_parity[i]) ) {
            (*num_known_SB) ++;
        } else {
            (*num_unk_und_SB) ++;
        }
    }
    /* count SC and allenes */
    for ( i = 0; i < Stereo->nNumberOfStereoCenters; i ++ ) {
        if ( !(nAtNumber = Stereo->nNumber[i]) || nAtNumber > pInChI->nNumberOfAtoms ) {
            return RI_ERR_PROGR; /* wrong data, should never happen */
        }
        if ( ATOM_PARITY_WELL_DEF(Stereo->t_parity[i]) ) {
            (*num_known_SC) ++;
        } else {
            (*num_unk_und_SC) ++;
        }
        el_number = pInChI->nAtom[nAtNumber-1];
        if ( el_number != el_number_P && el_number != el_number_As ) {
            continue;
        }
        ret = GetNumNeighborsFromInchi( pInChI, nAtNumber );
        if ( ret < 0 ) {
            return ret;
        }
        if ( 3 == ret ) {
            *num_SC_PIII  += (el_number_P  == el_number);
            *num_SC_AsIII += (el_number_As == el_number);
        }
    }
    return 2; /* Has Stereo */
}
/**************************************************************************************/
int bInpInchiComponentExists( InpInChI *pOneInput, int iInChI, int bMobileH, int k )
{
    if ( (INCHI_BAS != iInChI   && iInChI != INCHI_REC)  ||
         (TAUT_NON  != bMobileH && TAUT_YES != bMobileH) || k < 0 ) {
        return 0;
    }
    return ( k < pOneInput->nNumComponents[iInChI][bMobileH] &&
                 pOneInput->pInpInChI[iInChI][bMobileH] &&
                 pOneInput->pInpInChI[iInChI][bMobileH][k].nNumberOfAtoms > 0 &&
                !pOneInput->pInpInChI[iInChI][bMobileH][k].bDeleted );
}
/**************************************************************************************/
int bInpInchiComponentDeleted( InpInChI *pOneInput, int iInChI, int bMobileH, int k )
{
    if ( (INCHI_BAS != iInChI   && iInChI != INCHI_REC)  ||
         (TAUT_NON  != bMobileH && TAUT_YES != bMobileH) || k < 0 ) {
        return 0;
    }
    return ( k < pOneInput->nNumComponents[iInChI][bMobileH] &&
                 pOneInput->pInpInChI[iInChI][bMobileH] &&
                 pOneInput->pInpInChI[iInChI][bMobileH][k].nNumberOfAtoms > 0 &&
                 pOneInput->pInpInChI[iInChI][bMobileH][k].bDeleted );
}
/**************************************************************************************/
int bRevInchiComponentExists( StrFromINChI *pStruct, int iInChI, int bMobileH, int k )
{
    if ( !pStruct || /*!pStruct->at2 ||*/ !pStruct->num_atoms ||
         (INCHI_BAS != iInChI   && iInChI != INCHI_REC)  ||
         (TAUT_NON  != bMobileH && TAUT_YES != bMobileH) || k < 0 ) {
        return 0;
    }
    return ( k < pStruct->RevInChI.num_components[iInChI] &&
                 pStruct->RevInChI.pINChI[iInChI] &&
                 pStruct->RevInChI.pINChI[iInChI][k][bMobileH] &&
                 pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->nNumberOfAtoms > 0 &&
                !pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->bDeleted );
}

/**************************************************************************************/
int bRevInchiComponentDeleted( StrFromINChI *pStruct, int iInChI, int bMobileH, int k )
{
    if ( !pStruct || /*!pStruct->at2 ||*/ !pStruct->num_atoms ||
         (INCHI_BAS != iInChI   && iInChI != INCHI_REC)  ||
         (TAUT_NON  != bMobileH && TAUT_YES != bMobileH) || k < 0 ) {
        return 0;
    }
    return ( k < pStruct->RevInChI.num_components[iInChI] &&
                 pStruct->RevInChI.pINChI[iInChI] &&
                 pStruct->RevInChI.pINChI[iInChI][k][bMobileH] &&
                 pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->nNumberOfAtoms > 0 &&
                 pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->bDeleted );
}

/**************************************************************************************/
int DetectInpInchiCreationOptions ( InpInChI *pOneInput, int *bHasReconnected, int *bHasMetal,
                                             int *bHasFixedH, int *nModeFlagsStereo, int *bTautFlagsStereo )
{
    int ret=0, bHasStereo;
    int nModeFlagsValue=0, bTautFlagsValue; /* stereo flags */
    int iInChI, iMobileH, bIso, k, max_components, num_components;
    INChI *pInChI;
    int num_known_SB /*Stereo Bonds & Cumulenes >C==C==C==C< */;
    int num_known_SC /* Stereo Centers & Allenes >C=C=C< */;
    int num_unk_und_SB, num_unk_und_SC;
    int num_SC_PIII, num_SC_AsIII; /* has Phosphine or Arsine stereo center(s) */

    *bHasReconnected = *bHasFixedH = *nModeFlagsStereo = *bTautFlagsStereo = 0;
    nModeFlagsValue = bTautFlagsValue = bHasStereo = 0;
    num_known_SB = num_known_SC = num_unk_und_SB = num_unk_und_SC = num_SC_PIII = num_SC_AsIII = 0;
    *bHasMetal = 0;

    for ( iInChI   = 0; iInChI   < INCHI_NUM; iInChI ++   ) 
    {
        for ( iMobileH = 0; iMobileH < TAUT_NUM;  iMobileH ++ ) 
        {
            for ( bIso = 1; !nModeFlagsValue && 0 <= bIso; bIso -- ) 
            {
                switch( pOneInput->s[iInChI][iMobileH][bIso] ) 
                {
                case 1: /* SABS */
                    nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO;
                    break;
                case 2:
                    nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RELATIVE_STEREO;
                    break;
                case 3:
                    nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RACEMIC_STEREO;
                }
            }
        
            max_components = pOneInput->pInpInChI[iInChI][iMobileH]?
                             pOneInput->nNumComponents[iInChI][iMobileH] : 0;

            for ( k = num_components = 0; k < max_components; k ++ ) 
            {
                pInChI = pOneInput->pInpInChI[iInChI][iMobileH] + k;
                ret = CountStereoTypes(pInChI, 
                                       &num_known_SB, &num_known_SC,
                                       &num_unk_und_SB, &num_unk_und_SC,
                                       &num_SC_PIII, &num_SC_AsIII);
                if ( ret < 0 ) 
                {
                    return ret; /* error */
                }
                bHasStereo     += (ret == 2);
                if ( (ret > 0) ) 
                {
                    /* ret == 0 => Empty InChI, 1=> No Stereo, 2=> Has Stereo */
                    num_components ++;
                    *bHasReconnected |= ( iInChI   == INCHI_REC );
                    *bHasFixedH      |= ( iMobileH == TAUT_NON  );   
                }
                *bHasMetal |= bInChIHasReconnectedMetal( pInChI );
            }
        }
    }

    if ( (nModeFlagsValue & REQ_MODE_RELATIVE_STEREO) && (nModeFlagsValue & REQ_MODE_RACEMIC_STEREO) ) 
    {
        return RI_ERR_SYNTAX;
    }
    if ( bHasStereo && !nModeFlagsValue ) /* REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
    {  
        /* inversion does not change the stereo or no stereo at all */
        nModeFlagsValue = REQ_MODE_STEREO | REQ_MODE_ISO_STEREO; /* Abs */
    }

    if ( !num_known_SB && num_unk_und_SB ) 
    {
        ; /* full SUU option or SB part of it */
    } else 
    {
        nModeFlagsValue |= REQ_MODE_SB_IGN_ALL_UU; /* ignore Unknown/Undefind SB if no well-defined SB exist */
    }

    if ( !num_known_SC && num_unk_und_SC ) 
    {
        ; /* full SUU option or SB part of it */
    } else 
    {
        nModeFlagsValue |= REQ_MODE_SC_IGN_ALL_UU; /* ignore Unknown/Undefind SC if no well-defined SB exist */
    }

    /* Phosphine and Arsine Stereo */
    if ( num_SC_PIII ) 
    {
        bTautFlagsValue |= TG_FLAG_PHOSPHINE_STEREO;
    }
    /* Phosphine and Arsine Stereo */
    if ( num_SC_AsIII ) 
    {
        bTautFlagsValue |= TG_FLAG_ARSINE_STEREO;
    }

    *nModeFlagsStereo = nModeFlagsValue;
    *bTautFlagsStereo = bTautFlagsValue;
    
    return 0;
}

/******************************************************************************************************/
int bInChIHasReconnectedMetal( INChI *pInChI )
{
    int i;
    if ( pInChI && !pInChI->bDeleted && pInChI->nNumberOfAtoms && pInChI->nAtom ) 
    {
        for ( i = 0; i < pInChI->nNumberOfAtoms; i ++ ) 
        {
            if ( is_el_a_metal( (int)pInChI->nAtom[i] ) ) 
            {
                if ( pInChI->nNumberOfAtoms > 1 || (pInChI->nNum_H && pInChI->nNum_H[0]) )
                {
                    return 1;
                }
            }
        }
    }
    return 0;
}

/*****************************************************************************************/
int SetProtonsAndXchgIsoH( int bInChI2Structure, int bReqSplitOutputInChI, int bReqProtonsForEachComponent,
                           int bReqNonTaut, int bReqStereo, int num_components[INCHI_NUM],
                           MODE_PIXH nModeProtonIsoExchgH[INCHI_NUM], InpInChI *OneInput )
{
    int      j, k, k1, ret2=0, iINChI;
    int  bAvailableProtonsForEachComponent, bAvailableProtonsTotal;
    
    INCHI_HEAPCHK

    num_components[INCHI_BAS] = num_components[INCHI_REC] = 0;

    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) {
        nModeProtonIsoExchgH[iINChI] = MODE_PIXH_UNDEFINED;
        /* are totals of /p and/or /i/h available ? */
        bAvailableProtonsTotal            = 0 != OneInput->nNumProtons[iINChI][TAUT_YES].nNumRemovedProtons;
        for ( k1 = 0; k1 < NUM_H_ISOTOPES; k1 ++ ) {
            bAvailableProtonsTotal |= 0 !=OneInput->nNumProtons[iINChI][TAUT_YES].nNumRemovedIsotopicH[k1];
        }
        /* are /p and/or /i/h available for each component ? */
        bAvailableProtonsForEachComponent = (NULL != OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons);
        /* decision: add /p to each component, add total to the 1st, add total as one more component */
        /* In case of bInChI2Structure just keep totals if not available for each component */
        if ( bInChI2Structure ) {
            nModeProtonIsoExchgH[iINChI] = bAvailableProtonsForEachComponent? 
                                       MODE_PIXH_ADD_TO_EACH:
                                       MODE_PIXH_KEEP_TOTALS;
        } else
        if ( !bReqSplitOutputInChI ) {
            nModeProtonIsoExchgH[iINChI] = bAvailableProtonsForEachComponent?
                                       MODE_PIXH_ADD_TO_EACH :
                                       MODE_PIXH_ADD_TO_FIRST;
        } else
        if ( !bAvailableProtonsForEachComponent ) {
            nModeProtonIsoExchgH[iINChI] = bAvailableProtonsTotal?
                                       MODE_PIXH_ADD_A_PIXH_COMPONENT :
                                       MODE_PIXH_ADD_TO_FIRST;
        } else
        /* bAvailableProtonsForEachComponent && bReqSplitOutputInChI */
        if ( bReqProtonsForEachComponent ) {
            nModeProtonIsoExchgH[iINChI] =     MODE_PIXH_ADD_TO_EACH;
        } else {
            nModeProtonIsoExchgH[iINChI] = bReqNonTaut?
                                       MODE_PIXH_ADD_TO_EACH :
                                       MODE_PIXH_ADD_A_PIXH_COMPONENT;
        }
        /* remove unneeded data: protons for each component */
        if ( bAvailableProtonsForEachComponent &&
             nModeProtonIsoExchgH[iINChI] != MODE_PIXH_ADD_TO_EACH ) {
            inchi_free( OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons );
            OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons = NULL;
            bAvailableProtonsForEachComponent = 0;
        }
        /* remove unneeded data: total protons all components */
        if ( bAvailableProtonsTotal && nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_EACH ) {
            OneInput->nNumProtons[iINChI][TAUT_YES].nNumRemovedProtons = 0;
            for ( k1 = 0; k1 < NUM_H_ISOTOPES; k1 ++ ) {
                OneInput->nNumProtons[iINChI][TAUT_YES].nNumRemovedIsotopicH[k1] = 0;
            }
            bAvailableProtonsTotal = 0;
        }
        /* remove unneeded data: Fixed-H InChI; no protons data exist for Fixed-H */
        if ( !bReqNonTaut && OneInput->nNumComponents[iINChI][TAUT_NON] ) {
            j = TAUT_NON;
            for ( k = 0; k < OneInput->nNumComponents[iINChI][j]; k ++ ) {
                Free_INChI_Members( &OneInput->pInpInChI[iINChI][j][k] );
            }
            inchi_free( OneInput->pInpInChI[iINChI][j] );
            OneInput->pInpInChI[iINChI][j] = NULL;
            OneInput->nNumComponents[iINChI][j] = 0;
        }
#ifdef NEVER
        /* remove unneeded data: Mobile-H InChI ????? */
        if ( bReqNonTaut && OneInput->nNumComponents[iINChI][TAUT_NON] ) {
            j = TAUT_YES;
            for ( k = 0; k < OneInput->nNumComponents[iINChI][j]; k ++ ) {
                Free_INChI_Members( &OneInput->pInpInChI[iINChI][j][k] );
            }
            inchi_free( OneInput->pInpInChI[iINChI][j] );
            OneInput->pInpInChI[iINChI][j] = NULL;
            OneInput->nNumComponents[iINChI][j] = 0;
            nModeProtonIsoExchgH[iINChI] = MODE_PIXH_UNDEFINED;
            if ( OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons ) {
                inchi_free( OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons);
                OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons = NULL;
            }
        }
#endif
        /* add one more component containing only /p and /i/h */
        if ( (nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_A_PIXH_COMPONENT &&
             OneInput->nNumComponents[iINChI][TAUT_YES])       ||
             /* always add one deleted component if no non-taut InChI is available */
             (bInChI2Structure && !bAvailableProtonsForEachComponent &&
             !OneInput->nNumComponents[iINChI][TAUT_NON] &&
             OneInput->nNumComponents[iINChI][TAUT_YES]) ) {
            int nPrevLen, nLen=0;
            j = TAUT_YES;
            nPrevLen = OneInput->nNumComponents[iINChI][j];
            for ( k = 0; k < nPrevLen; k ++ ) {
                nLen += !OneInput->pInpInChI[iINChI][j][k].bDeleted;    
            }
            if ( nLen == nPrevLen ) {
                /* add one more component */
                INChI *pInChI = (INChI *)inchi_calloc( nLen+1, sizeof(*pInChI) );
                if ( !pInChI ) {
                    ret2 = RI_ERR_ALLOC;
                    goto exit_error;
                }
                memcpy( pInChI, OneInput->pInpInChI[iINChI][j], nLen*sizeof(*pInChI) );
                inchi_free( OneInput->pInpInChI[iINChI][j] );
                OneInput->pInpInChI[iINChI][j] = pInChI;
            }
            OneInput->nNumComponents[iINChI][j] = nLen+1;
            
            for ( k = nLen; k < nPrevLen; k ++ ) {
                Free_INChI_Members( &OneInput->pInpInChI[iINChI][j][k] );
                memset( &OneInput->pInpInChI[iINChI][j][k], 0, sizeof(OneInput->pInpInChI[iINChI][j][k]));
            }
            /* mark the last component as a proton */
            if ( 0 > (ret2 = nFillOutProtonMobileH( OneInput->pInpInChI[iINChI][j]+nLen ) ) ) {
                goto exit_error;
            }
        }
        INCHI_HEAPCHK
       
        /* remove unneeded Stereo and/or Fixed H */
        if ( !bReqStereo ) {
            for ( j = 0; j < TAUT_NUM; j ++ ) {
                for ( k = 0; k < OneInput->nNumComponents[iINChI][j]; k ++ ) {
                    if ( OneInput->pInpInChI[iINChI][j][k].Stereo ) {
                        Free_INChI_Stereo(OneInput->pInpInChI[iINChI][j][k].Stereo);
                        inchi_free( OneInput->pInpInChI[iINChI][j][k].Stereo );
                        OneInput->pInpInChI[iINChI][j][k].Stereo = NULL;
                    }
                    if ( OneInput->pInpInChI[iINChI][j][k].StereoIsotopic ) {
                        Free_INChI_Stereo(OneInput->pInpInChI[iINChI][j][k].StereoIsotopic);
                        inchi_free( OneInput->pInpInChI[iINChI][j][k].StereoIsotopic );
                        OneInput->pInpInChI[iINChI][j][k].StereoIsotopic = NULL;
                    }
                    INCHI_HEAPCHK
                }
            }
        }
    }

    num_components[INCHI_BAS] = inchi_max( OneInput->nNumComponents[INCHI_BAS][TAUT_YES],
                                           OneInput->nNumComponents[INCHI_BAS][TAUT_NON] );
    num_components[INCHI_REC] = inchi_max( OneInput->nNumComponents[INCHI_REC][TAUT_YES],
                                           OneInput->nNumComponents[INCHI_REC][TAUT_NON] );

exit_error:
    return ret2;
}
/******************************************************************************************/
int GetInChIFormulaNumH( INChI *pInChI, int *nNumH )
{  /* get number of H including bridging hydrogen atoms */
    const char *p, *q;
    *nNumH = 0;
    if ( pInChI->szHillFormula ) {
        for ( p = strchr( pInChI->szHillFormula, 'H'); p; p = strchr(p, 'H') ) {
            p ++;
            if ( !islower( UCINT *p ) ) {
                /* found hydrogen in the formula */
                if ( isdigit( UCINT *p ) ) {
                    *nNumH += (int)inchi_strtol( p, &q, 10 );
                    p = q;
                } else {
                    *nNumH += 1;
                }
            }
        }
    }
    return 0;
}
/******************************************************************************************/
int GetInChINumH( INChI *pInChI, int *nNumH )
{
    int i, j, nNumTautGroups, iTautGroup, nTautGroupLen, lenTautomer;
    *nNumH = 0;
    for ( i = 0; i < pInChI->nNumberOfAtoms; i ++ ) {
        *nNumH += ( pInChI->nAtom[i] == EL_NUMBER_H ); /* bridging H */
        *nNumH += pInChI->nNum_H[i];
    }
    /* earlier nNum_H_fixed[] should have been added to pInChI->nNum_H[] */
    /*
    if ( pInChI->nNum_H_fixed ) {
        for ( i = 0; i < pInChI->nNumberOfAtoms; i ++ ) {
            *nNumH += pInChI->nNum_H_fixed[i];
        }
    }
    */
    if ( pInChI->lenTautomer > 3 && pInChI->nTautomer ) {
        lenTautomer = pInChI->lenTautomer;
        j = 0;
        nNumTautGroups = pInChI->nTautomer[j ++];
        for ( iTautGroup = 0; j < lenTautomer && iTautGroup < nNumTautGroups; iTautGroup ++, j += nTautGroupLen ) {
            nTautGroupLen = pInChI->nTautomer[j]+1;
            *nNumH += pInChI->nTautomer[j+1];
        }
        if ( iTautGroup != nNumTautGroups || j != lenTautomer ) {
            return RI_ERR_PROGR;
        }
    }
    if ( pInChI->nNum_H_fixed && (pInChI->lenTautomer || pInChI->nTautomer) ) {
        return RI_ERR_PROGR;
    }
    return 0;
}
/******************************************************************************************/
int GetInChIIsoH( INChI *pInChI, int nNumIsotopicH[NUM_H_ISOTOPES] )
{
    int i;
    for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
        nNumIsotopicH[i] = 0;
    }
    for ( i = 0; i < pInChI->nNumberOfIsotopicAtoms; i ++ ) {
        if ( pInChI->IsotopicAtom[i].nIsoDifference > 0 &&
             pInChI->IsotopicAtom[i].nIsoDifference <= NUM_H_ISOTOPES ) {
            if ( !pInChI->nAtom ||
                 !pInChI->IsotopicAtom[i].nAtomNumber ||
                 pInChI->IsotopicAtom[i].nAtomNumber > pInChI->nNumberOfAtoms ) {
                return RI_ERR_PROGR;
            }
            if ( pInChI->nAtom[pInChI->IsotopicAtom[i].nAtomNumber-1] == EL_NUMBER_H ) {
                /* isotopic H in connection table */
                nNumIsotopicH[ pInChI->IsotopicAtom[i].nIsoDifference-1 ] ++;
            }
        }
        nNumIsotopicH[0] += pInChI->IsotopicAtom[i].nNum_H;
        nNumIsotopicH[1] += pInChI->IsotopicAtom[i].nNum_D;
        nNumIsotopicH[2] += pInChI->IsotopicAtom[i].nNum_T;
    }
    return 0;
}
/******************************************************************************************/
typedef struct tagNumElem 
{
    int num;
    /*
    int iso;
    */
} NUM_ELEM;


/***************************************************************************************/
int InChILine2Data(INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, 
                   char **pStr, int *pState, int *nErr,
                   INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                   int nNumComponents[INCHI_NUM][TAUT_NUM],
                   REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                   int s[INCHI_NUM][TAUT_NUM][2], 
                   int bReadCoord, int bInchi2Struct, INCHI_MODE nMode,
                   int *bStdFormat, int *bInputHasSaveOpt, unsigned char *inp_save_opt_bits)
{
    int iINChI, i, j, k, m, len1, len2, ret2 = 0, retAux = 0, stateAux = 0;
    int ret, tot_charge[INCHI_NUM][TAUT_NUM];
    int i1, i2, i3;
    int kc;
    NUM_ELEM *num_elem[INCHI_NUM][TAUT_NUM];


#if ( FIX_I2I_STEREOCONVERSION_BUG == 1 )
/* (2008-03-06)   1=> Fix bug of i2i conversion SAbs-->(SRel||Srac) */
/*                    (converter does not placed proper stereo to output) */
    
    /* set new stereo type as requested by conversion option */
    int target_stereo_type = 1;                
    if (nMode & REQ_MODE_RELATIVE_STEREO)
        target_stereo_type = 2;
    else if (nMode & REQ_MODE_RACEMIC_STEREO)
        target_stereo_type = 3;
#endif


    memset( num_elem, 0, sizeof(num_elem) );

    ret = ReadInChILine(pInp, pLine, pStr, pState, pInpInChI, 
                        nNumComponents, nNumProtons, s,
                        bStdFormat, bInputHasSaveOpt, inp_save_opt_bits);


#if ( FIX_I2I_STEREOCONVERSION_BUG == 1 )
    /* modify stereo type for layers as requested */
    if (target_stereo_type > 1)
        for (i1=0;i1<INCHI_NUM;i1++)
            for (i2=0;i2<TAUT_NUM;i2++)
                for (i3=0;i3<2;i3++)
                    if (s[i1][i2][i3] != 0) 
                    {
                        if ( target_stereo_type!=1) 
                            /* do not allow conversion SRel=>SAbs, SRac=>SAbs */
                            s[i1][i2][i3] = target_stereo_type;						
                    }
#endif


    *nErr = 0;
    if ( (ret == RI_ERR_EOL) &&
                       nNumComponents[INCHI_BAS][TAUT_YES] 
                     + nNumComponents[INCHI_BAS][TAUT_NON] && bReadCoord ) {
        retAux = ReadInChICoord( pInp, pLine, &stateAux, pInpInChI, nNumComponents );
    }
    if ( (ret == RI_ERR_EOL || ret == RI_ERR_EOF) &&
                       nNumComponents[INCHI_BAS][TAUT_YES] 
                     + nNumComponents[INCHI_BAS][TAUT_NON] ) {
        /* post-processing: add omitted layers */
        *pState = IST_MATERIAL_BALANCE_ERROR;
        for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) 
        {
            
            for ( j = 0; j < TAUT_NUM; j ++ ) 
            {
                /* for Mobile/Fixed H (j) ... */

                int bIsotopic, bStereoType, bStereoTypeAlt;
                int nMH2FH_AltInv=0, nFH2iFH_AltInv=0 /*, niMH2iFH_AltInv=0, nMH2iMH_AltInv=0*/;
                int jAlt = ALT_TAUT(j);
                INCHI_MODE  nFlags = 0, nFlagsAlt = 0;
                /* get stereo type: ABS, REL, RAC, or nothing */
                tot_charge[iINChI][j] = 0;
                for ( bIsotopic = bStereoType = bStereoTypeAlt = 0; bIsotopic < 2; bIsotopic ++ ) {
                    if ( !bStereoType || bStereoType < s[iINChI][j][bIsotopic] ) {
                        bStereoType = s[iINChI][j][bIsotopic];
                    }
                    if ( !bStereoTypeAlt || bStereoTypeAlt < s[iINChI][jAlt][bIsotopic] ) {
                        bStereoTypeAlt = s[iINChI][jAlt][bIsotopic];
                    }
                    nFlags = bStereoType      ==2? INCHI_FLAG_REL_STEREO : bStereoType   ==3? INCHI_FLAG_RAC_STEREO : 0;
                    nFlagsAlt = bStereoTypeAlt==2? INCHI_FLAG_REL_STEREO : bStereoTypeAlt==3? INCHI_FLAG_RAC_STEREO : 0;
                }
                /* set stereo type to each component */
                /* add missing nNum_H and nConnTable */
                if ( nNumComponents[iINChI][j] ) {
                    num_elem[iINChI][j] = (NUM_ELEM *)inchi_calloc( nElDataLen+1, sizeof(num_elem[0][0][0]) );
                    if ( !num_elem[iINChI][j] ) {
                        ret2 = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                for ( k = 0; k < nNumComponents[iINChI][j]; k ++ ) 
                {
                    /* for each component k ... */

                    if ( pInpInChI[iINChI][j] ) {
                        INChI *pInChI = &pInpInChI[iINChI][j][k];
                        INChI *pInChI_Alt = (k<nNumComponents[iINChI][jAlt] &&
                                             pInpInChI[iINChI][jAlt] &&
                                             /*pInpInChI[iINChI][jAlt]->nNumberOfAtoms)? pInpInChI[iINChI][jAlt]:NULL;*/ /* 2007-09-25 DT */
                                             pInpInChI[iINChI][jAlt][k].nNumberOfAtoms)? &pInpInChI[iINChI][jAlt][k]:NULL;
                                             
                        if ( nFlags ) {
                            pInChI->nFlags |= nFlags;
                        } else
                        if ( j == TAUT_NON && !nFlags && nFlagsAlt ) {
                            pInChI->nFlags |= nFlagsAlt;
                        }
                        /**** add empty immobile H (nNum_H) if it is missing ****/
                        if ( !pInChI->nNum_H &&
                             !(pInChI->nNum_H = (S_CHAR *)inchi_calloc( pInChI->nNumberOfAtoms+1, sizeof(pInChI->nNum_H[0]) ) ) ) {
                            ret2 = RI_ERR_ALLOC;
                            goto exit_function;
                        }
                        /**** add single atom nConnTable if it is missing ****/
                        if ( !pInChI->nConnTable ) {
                            AT_NUMB *pCT;
                            int      lenCT;
                            if ( j == TAUT_NON && k < nNumComponents[iINChI][TAUT_YES] &&
                                 (pCT = pInpInChI[iINChI][TAUT_YES][k].nConnTable) &&
                                 (lenCT = pInpInChI[iINChI][TAUT_YES][k].lenConnTable) > 0) {
                                if ( !(pInChI->nConnTable = (AT_NUMB *)inchi_calloc( lenCT+1, sizeof(pInChI->nConnTable[0]) ) ) ) {
                                    ret2 = RI_ERR_ALLOC;
                                    goto exit_function;
                                }
                                memcpy( pInChI->nConnTable, pCT, lenCT*sizeof(pInChI->nConnTable[0]) );
                                pInChI->lenConnTable = lenCT;
                            } else {
                                if ( j == TAUT_YES && pInChI->nNumberOfAtoms > 1 ) {
                                    *pState = IST_MOBILE_H_CONNECTIONS + (iINChI==INCHI_REC? IST_HAPPENED_IN_RECMET : 0);
                                    ret2 = RI_ERR_SYNTAX;
                                    goto exit_function;
                                }
                                if ( !(pInChI->nConnTable = (AT_NUMB *)inchi_calloc( pInChI->nNumberOfAtoms+1, sizeof(pInChI->nConnTable[0]) ) ) ) {
                                    ret2 = RI_ERR_ALLOC;
                                    goto exit_function;
                                }
                                pInChI->lenConnTable = 1;
                                pInChI->nConnTable[0] = 1;
                            }
                        } else
                        if ( pInChI->nConnTable && !pInChI->lenConnTable && pInChI->nNumberOfAtoms == 1 ) {
                            pInChI->nConnTable[0] = 1;
                            pInChI->lenConnTable  = 1;
                        }
                        /**** copy charge: Mobile H --> Fixed H; ****/ 
                        if ( j == TAUT_NON ) {
                            /*
                            if ( pInChI->nTotalCharge == NO_VALUE_INT ) {
                                pInChI->nTotalCharge = 0;
                            } else
                            */
                            if ( !pInChI->nTotalCharge && k < nNumComponents[iINChI][TAUT_YES] ) {
                                INChI *pAltInChI = &pInpInChI[iINChI][TAUT_YES][k]; /* Mobile H InChI */
                                if ( pAltInChI->nTotalCharge && pAltInChI->nTotalCharge != NO_VALUE_INT ) {
                                    pInChI->nTotalCharge = pAltInChI->nTotalCharge;
                                }
                            }
                        }
                        /***** Fixed H: add pInChI->nNum_H_fixed to pInChI->nNum_H ****/
                        if ( j == TAUT_NON && pInChI->nNum_H && pInChI->nNum_H_fixed ) {
                            for ( m = 0; m < pInChI->nNumberOfAtoms; m ++ ) {
                                pInChI->nNum_H[m] += pInChI->nNum_H_fixed[m];
                            }
                        }
                        /***** copy isotopic atoms: Mobile H --> Fixed H ******/ 
                        if ( j == TAUT_YES && pInChI->nNumberOfIsotopicAtoms &&
                             k < nNumComponents[iINChI][TAUT_NON] ) {
                            INChI *pAltInChI = &pInpInChI[iINChI][TAUT_NON][k]; /* Fixed H InChI */

                            if ( !pAltInChI->nNumberOfIsotopicAtoms ) {
                                ret2=CopySegment( pAltInChI, pInChI, CPY_ISO_AT, 0, 0);
                                if ( ret2 < 0 ) {
                                    goto exit_function;
                                }
                            }
                        }
                        /**** copy coordinates: Mobile H --> Fixed H ******/
                        if ( j == TAUT_YES && pInChI->IsotopicTGroup &&
                             k < nNumComponents[iINChI][TAUT_NON] ) {
                            INChI *pAltInChI = &pInpInChI[iINChI][TAUT_NON][k]; /* Fixed H InChI */

                            if ( !pAltInChI->IsotopicTGroup ) {
                                XYZ_COORD *pxyz = (XYZ_COORD *)inchi_calloc( pInChI->nNumberOfAtoms, sizeof(pxyz[0]));
                                if ( pxyz ) {
                                    memcpy( pxyz, pInChI->IsotopicTGroup, pInChI->nNumberOfAtoms * sizeof(pxyz[0]) );
                                    pAltInChI->IsotopicTGroup = (INChI_IsotopicTGroup *)pxyz;
                                } else {
                                    ret2 = RI_ERR_ALLOC;
                                    goto exit_function;
                                }
                            }
                        }

                        /********************************************************
                         *                                                      *
                         *            Restore omitted stereo seqments           *
                         *                                                      *
                         * order of restoring:                                  *
                         *                                                      *
                         * 1. Fixed H            (F) -> (FI) Isotopic Fixed H   *
                         * 2. Mobile H           (M) -> (F)  Fixed H            *
                         * 3. Isotopic Mobile H (MI) -> (FI) Isotopic Fixed H   *
                         * 4. Mobile H           (M) -> (MI) Isotopic Mobile H  *
                         *                                                      *
                         ********************************************************/

                        /***** (4) copy stereo: Mobile H --> isotopic Mobile H ******/
                        if ( j == TAUT_YES ) {
                            int bIso = pInChI->nNumberOfIsotopicAtoms || 
                                   (pInChI->StereoIsotopic &&
                                    pInChI->StereoIsotopic->nNumberOfStereoCenters 
                                  + pInChI->StereoIsotopic->nNumberOfStereoBonds) ||
                                  (pInChI_Alt && pInChI_Alt->nNumberOfIsotopicAtoms);

                            /* non-isotopic Mobile H => isotopic Mobile H */
                            if ( bIso ) {
                                if ( pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                     (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->t_parity) ) {
                                    if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3, 1, 0)) ||
                                         ((!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) &&
                                         0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3_M, 1, 0)))) {
                                        goto exit_function;
                                    }
                                    if ( (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                         if ( pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT ) {
                                             pInChI->Stereo->nCompInv2Abs = s[iINChI][j][0]>0? 2 : 0;
                                         }
                                         if ( pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                             pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][j][1]>0? 2 : 0;
                                         }
                                    }
                                } else
                                /* copy sp3 inversion info: non-isotopic Mobile H => isotopic Mobile H  */
                                if ( pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                     pInChI->StereoIsotopic && pInChI->StereoIsotopic->nNumberOfStereoCenters &&
                                     pInChI->Stereo->nCompInv2Abs ) {
                                    if ( (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
                                         pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT &&
                                         pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                        pInChI->Stereo->nCompInv2Abs = s[iINChI][j][0]>0? 2 : 0;
                                        pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][j][1]>0? 2 : 0;
                                    } else
                                    if (!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs ) {
                                        pInChI->StereoIsotopic->nCompInv2Abs = pInChI->Stereo->nCompInv2Abs;
                                    }
                                }
                            }
                            if ( bIso &&
                                 pInChI->Stereo && pInChI->Stereo->nNumberOfStereoBonds &&
                                 (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->b_parity) ) {
                                if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP2, 1, 0)) ) {
                                    goto exit_function;
                                }
                            }
                        }
                        /***** (0) set nCompInv2Abs to Fixed-H *********************************/
                        if ( j == TAUT_NON ) {
                            if ( pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                 pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT ) {
                                /* case of /sN and /t... in non-isotopic Mobile-H, no /s in non-isotopic Fixed-H */
                                if ( !s[iINChI][j][0] && s[iINChI][jAlt][0]>0 &&  /* /sN is not present in F and is present in M */
                                     pInChI_Alt && pInChI_Alt->Stereo && pInChI_Alt->Stereo->nNumberOfStereoCenters ) {
                                    /* inherit from Mobile-H */
                                    /* /s1 in M and MI; /m1 or /m0 in MI; /m. in M; no /m in F. Inherit MI->FI. Added 10-15-2007 */
                                    if ( pInChI_Alt->Stereo->nCompInv2Abs == 0 &&                    /*  M: /m. ; means no /m for this component */
                                         pInChI->Stereo->nCompInv2Abs     == NO_VALUE_INT &&         /*  F: no /m segment for all components */
                                         pInChI_Alt->StereoIsotopic &&                               /*  MI: present */
                                         pInChI_Alt->StereoIsotopic->nCompInv2Abs != 0 &&
                                         pInChI_Alt->StereoIsotopic->nCompInv2Abs != NO_VALUE_INT && /* MI:    /m0 or /m1  */
                                         !s[iINChI][j][0] && !s[iINChI][j][1] &&                     /* F, FI: no /s       */
                                         s[iINChI][jAlt][0] == 1 && s[iINChI][jAlt][1] == 1          /* M, MI: /s1 and /s1 */
                                       ) {
                                        /* copy /m from MI to FI */
                                        if ( 0 > (ret2 = CopySegment( pInChI, pInChI_Alt, CPY_SP3_M, 1, 1)) ) {
                                            goto exit_function;
                                        }
                                    } else
                                    /* the following if(){...} was added to fix m1 bug 2007-09-25 DT */
                                    if ( pInChI_Alt->Stereo->nCompInv2Abs != NO_VALUE_INT && s[iINChI][jAlt][0] == 1 ) {
                                        pInChI->Stereo->nCompInv2Abs = pInChI_Alt->Stereo->nCompInv2Abs;
                                    } else
                                    /* M and MI contain /sN and /sN, N=2,3. Added 10-15-2007 */
                                    if ( pInChI_Alt->Stereo->nCompInv2Abs == NO_VALUE_INT &&
                                         pInChI->Stereo->nCompInv2Abs     == NO_VALUE_INT &&
                                         !s[iINChI][j][0] && !s[iINChI][j][1] &&
                                         (s[iINChI][jAlt][0] & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
                                         (s[iINChI][jAlt][1] & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                        int bIso = pInChI->nNumberOfIsotopicAtoms || 
                                               (pInChI->StereoIsotopic &&
                                                pInChI->StereoIsotopic->nNumberOfStereoCenters 
                                              + pInChI->StereoIsotopic->nNumberOfStereoBonds) ||
                                               (pInChI_Alt && pInChI_Alt->nNumberOfIsotopicAtoms);
                                        if ( bIso ){
                                            if ( !pInChI_Alt->StereoIsotopic &&  /* create zero/NULL-initialized pInChI_Alt->StereoIsotopic */
                                                 0 > (ret2 = CopySegment( pInChI_Alt, pInChI_Alt, CPY_SP3_M, 1, -1))) {
                                                goto exit_function;
                                            }
                                            pInChI_Alt->StereoIsotopic->nCompInv2Abs = 2;  /* MI: /m1 or /m0 */
                                            pInChI_Alt->Stereo->nCompInv2Abs         = 0;  /* M:  /m. ; no /m for this component */
                                            pInChI->Stereo->nCompInv2Abs = NO_VALUE_INT+1; /* FI: Stereo->CompInv2Abs=0, StereoIsotopic->CompInv2Abs=1 or -1 */
                                        } else {
                                            pInChI->Stereo->nCompInv2Abs     = 2; /* F:  /m1 or /m0, omitted from InChI as a repetition */
                                            pInChI_Alt->Stereo->nCompInv2Abs = 2; /* M:  /m1 or /m0; in Srel/SRac case the value = 2 */
                                        }
                                    } else {
                                        pInChI->Stereo->nCompInv2Abs = 2;         /* F:  /m1 or /m0, omitted from InChI as a repetition */
                                        pInChI_Alt->Stereo->nCompInv2Abs = 2;     /* M:  /m1 or /m0; in Srel/SRac case the value = 2 */
                                    }
                                } else
                                /* case of /sN in Isotopic Fixed-H only, /t... in Fixed-H, no /m (2007-08-27 DT) */
                                if ( !s[iINChI][j][0] && !s[iINChI][jAlt][0] && /* /sN in Fixed-H isotopic only */
                                     (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
                                     !(pInChI->StereoIsotopic && pInChI->StereoIsotopic->nNumberOfStereoCenters) &&
                                     /*!(pInChI_Alt && pInChI_Alt->Stereo && pInChI_Alt->Stereo->nNumberOfStereoCenters) &&*/
                                     !(pInChI_Alt && pInChI_Alt->StereoIsotopic && pInChI_Alt->StereoIsotopic->nNumberOfStereoCenters) ) {
                                    pInChI->Stereo->nCompInv2Abs = NO_VALUE_INT+1; /* Stereo->CompInv2Abs=0, StereoIsotopic->CompInv2Abs=1 or -1 */
                                } else {
                                    pInChI->Stereo->nCompInv2Abs = s[iINChI][j][0]>0? 2 : 0;
                                }
                            }
                        }
                        /***** (1) copy stereo: non-isotopic Fixed H --> isotopic Fixed H ******/
                        if ( j == TAUT_NON ) {
                            int bIso = pInChI->nNumberOfIsotopicAtoms || 
                                   (pInChI->StereoIsotopic &&
                                    pInChI->StereoIsotopic->nNumberOfStereoCenters 
                                  + pInChI->StereoIsotopic->nNumberOfStereoBonds) ||
                                   (pInChI_Alt && pInChI_Alt->nNumberOfIsotopicAtoms);
                            /* non-isotopic Fixed H => isotopic Fixed H */
                            if ( bIso ) {
                                if ( pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                    (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->t_parity) ) {
                                    /* -- replaced 2007-08-27 by (aaa), see below -- DT
                                    if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3, 1, 0)) ||
                                         !(pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) &&
                                         0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3_M, 1, 0))) {
                                        goto exit_function;
                                    }
                                    */
                                    /*----------- replacement (aaa) begin 2007-08-27 DT */
                                    if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3, 1, 0)) ) {
                                        goto exit_function;
                                    }
                                    if ( pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT+1 ) {
                                        pInChI->Stereo->nCompInv2Abs = 0;
                                        pInChI->StereoIsotopic->nCompInv2Abs = 2;
                                    } else
                                    if ( !(pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) &&
                                         0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3_M, 1, 0))) {
                                        goto exit_function;
                                    }
                                    /*----------- replacement (aaa) end 2007-08-27 DT */
                                    if ( (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                         if ( pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT ) {
                                             pInChI->Stereo->nCompInv2Abs = s[iINChI][j][0]>0? 2 : 0;
                                         }
                                         if ( pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                             pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][j][1]>0? 2 : 0;
                                         }
                                    }
#ifdef NEVER
                                    if ( (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
                                         !s[iINChI][j][0] && s[iINChI][j][0]>0 ) {
                                        /* copied Rel/Rac stereo to Iso; /s is in Iso /s is not in non-Iso */
                                        /* this means all difference in stereo is in inversion */
                                        if ( pInChI->Stereo->nCompInv2Abs         == NO_VALUE_INT &&
                                             pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                            pInChI->Stereo->nCompInv2Abs = 0;         /* missing */
                                            pInChI->StereoIsotopic->nCompInv2Abs = 2; /* unusual value */
                                        }
                                    }
#endif
                                } else
                                /* copy sp3 inversion info: non-isotopic Fixed H --> isotopic Fixed H  */
                                if ( pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                     pInChI->StereoIsotopic && pInChI->StereoIsotopic->nNumberOfStereoCenters &&
                                     pInChI->Stereo->nCompInv2Abs ) {
                                    if ( (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
                                         pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT &&
                                         pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                        pInChI->Stereo->nCompInv2Abs = s[iINChI][j][0]>0? 2 : 0;
                                        pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][j][1]>0? 2 : 0;
                                    } else
                                    if (!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) {
                                        pInChI->StereoIsotopic->nCompInv2Abs = pInChI->Stereo->nCompInv2Abs;
                                    }
                                }   
                            }
                            if ( bIso &&
                                 pInChI->Stereo && pInChI->Stereo->nNumberOfStereoBonds &&
                                 (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->b_parity) ) {
                                if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP2, 1, 0)) ) {
                                    goto exit_function;
                                }
                            }
                        }

                        /***** copy stereo: Mobile H --> Fixed H ******/
                        if ( j == TAUT_NON && k < nNumComponents[iINChI][TAUT_YES] ) {
                            INChI *pAltInChI = &pInpInChI[iINChI][TAUT_YES][k]; /* Mobile H InChI */
                            int bIso = pInChI->nNumberOfIsotopicAtoms || 
                                   (pInChI->StereoIsotopic &&
                                    pInChI->StereoIsotopic->nNumberOfStereoCenters 
                                  + pInChI->StereoIsotopic->nNumberOfStereoBonds) ||
                                   (pAltInChI && (
                                   pAltInChI->nNumberOfIsotopicAtoms || 
                                   (pAltInChI->StereoIsotopic &&
                                    pAltInChI->StereoIsotopic->nNumberOfStereoCenters 
                                  + pAltInChI->StereoIsotopic->nNumberOfStereoBonds)) );
                            int bNo_InChI_t = (!pInChI->Stereo || !pInChI->Stereo->t_parity);
                            int bNo_InChI_m = (!pInChI->Stereo || NO_VALUE_INT == pInChI->Stereo->nCompInv2Abs);

                            /* (2) non-isotopic Mobile H => non-isotopic Fixed H */
                            if ( pAltInChI->Stereo && pAltInChI->Stereo->nNumberOfStereoCenters &&
                                 (!pInChI->Stereo || !pInChI->Stereo->t_parity) ) 
                            {
#if ( FIX_I2I_STEREOCONVERSION_BUG2 == 1 )	
                                /* (2008-04-02)   1=> Fix bug of i2i conversion SAbs-->(SRel||Srac) */
                                /*                    (converter skipped empty '/t' or sometimes produced an excess one */

                                /* check whether t stereo is actually present */
                                int bHave_t_stereo = 1;
                                if (pInChI->Stereo)
                                    bHave_t_stereo = pInChI	->Stereo->nNumberOfStereoCenters;
                                /* account for stereobonds present */
                                if ( bHave_t_stereo < 1 )
                                    if ( pInChI->Stereo->nNumberOfStereoBonds > 0 )
                                        bHave_t_stereo=1;
                                /* copy stereo anyway ... */
#endif
                                if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3, 0, 0)) ||
                                     ((!pInChI->Stereo->nCompInv2Abs || NO_VALUE_INT == pInChI->Stereo->nCompInv2Abs) &&
                                     0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3_M, 0, 0))) )
                                {
                                    goto exit_function;
                                }

#if ( FIX_I2I_STEREOCONVERSION_BUG2 == 1 )				
                                /* ... correct just copied stereo if applicable */
                                if ( (s[iINChI][j][0] < 1)	&&
                                     (bHave_t_stereo <1 )		&&
                                     (pAltInChI->Stereo->nNumberOfStereoCenters > 0) &&
                                     (s[iINChI][jAlt][0] < 1) )
                                {
                                    /* (2010-02-28) if not all stereo centers are unknown/undefined */
                                    /*  at which condition stereo still should present .. */
                                    int all_UU = 1;
                                    for (kc=0; kc<pAltInChI->Stereo->nNumberOfStereoCenters; kc++)
                                    {
                                        if ( (pAltInChI->Stereo->t_parity[kc] != AB_PARITY_UNKN) &&
                                             (pAltInChI->Stereo->t_parity[kc] != AB_PARITY_UNDF) )
                                        {
                                            all_UU=0;
                                            break;
                                        }
                                    }
                                    if (!all_UU)
                                        pInChI->Stereo->nNumberOfStereoCenters = 0;
                                }					
#endif

                                /* in case of missing nCompInv2Abs, 2005-05-10 */
                                if ( (pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT) &&
                                     (nFlagsAlt & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                    if ( s[iINChI][jAlt][0] > 0 && s[iINChI][j][0] > 0 ) {
                                         /* suppose once in a while only non-taut stereo changes if inverted */
                                        pAltInChI->Stereo->nCompInv2Abs = (++nMH2FH_AltInv)%2? 2:0;
                                        pInChI->Stereo->nCompInv2Abs = 2;
                                    } else
                                    /* Mobile-H: /t.. /sN; Mobile-H isotopic: /sN (n=2 or 3), not /t...; Fixed-H layer is present, has no /t, no /i/t */
                                    /* Mobile-H /sN was caused by another component that would have same /mN in all layers */
                                    /* therefore, in case of Abs. Stereo, Mobile-H stereo isotopic stereo would have /m1 */
                                    /* In case of Rel/Rac stereo, since no /m1 could occur in Mobile-H isotopic, */
                                    /* no pAltInChI->StereoIsotopic or pInChI->StereoIsotopic have been created yet. */
                                    /* added 10-11-2007 to fix i2i bug for Rel/Rac stereo */
                                    if ( nNumComponents[iINChI][j] > 1 &&
                                         bNo_InChI_t && bNo_InChI_m /* no /t... or /mN in Fixed-H  */ && !nFlags &&
                                         !(pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->t_parity) &&
                                         !(pInChI->StereoIsotopic    && pInChI->StereoIsotopic->t_parity) &&
                                          s[iINChI][j][0]==0 && s[iINChI][j][1] == 0 &&
                                         /* /sN, N=2 or 3 only in Mobile-H AND Mobile-H isotopic */
                                         (s[iINChI][jAlt][0] & ((INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO))) &&
                                         (s[iINChI][jAlt][1] & ((INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO))) )
                                    {
                                        if ( bIso ) {
                                            /* create two zero/NULL-initialized isotopic stereo if they do not exist */
                                            if ( (!pInChI->StereoIsotopic    && 0 > (ret2 = CopySegment( pInChI,    pAltInChI, CPY_SP3_M, 1, -1)))
                                                /* -- the following will be created later, in TAUT_YES part of the code -- */
                                                || (!pAltInChI->StereoIsotopic && 0 > (ret2 = CopySegment( pAltInChI, pAltInChI, CPY_SP3_M, 1, -1))) ) {
                                                goto exit_function;
                                            }
                                            /* same value = 2 for MI and FI; here we assign only FI */
                                            pInChI->StereoIsotopic->nCompInv2Abs = 2;
                                            pInChI->Stereo->nCompInv2Abs = 0;
                                            /* -- the following will NOT be assigned later, in TAUT_YES part of the code -- */
                                            pAltInChI->StereoIsotopic->nCompInv2Abs = 2;
                                            pAltInChI->Stereo->nCompInv2Abs = 0;
                                            /* */
                                        } else {
                                            if ( NO_VALUE_INT == pInChI->Stereo->nCompInv2Abs &&
                                                 NO_VALUE_INT == pAltInChI->Stereo->nCompInv2Abs ) {
                                                pInChI->Stereo->nCompInv2Abs = 2;
                                                pAltInChI->Stereo->nCompInv2Abs = 2;
                                            }
                                        }
                                    } else
                                    if ( (s[iINChI][jAlt][0] > 0 || s[iINChI][j][0] > 0) && s[iINChI][j][0] >= 0 )
                                        pInChI->Stereo->nCompInv2Abs =  2;
                                    else
                                    /* Mobile-H: /t..., no /sN; Mobile-H isotopic: /s2 or /s3, not /t; Fixed-H layer is present, has no /t, no /i/t */
                                    /* therefore, in case of Abs. Stereo, Mobile-H stereo isotopic stereo would have /m1 */
                                    /* In case of Rel/Rac stereo, since no /m1 could occur in Mobile-H isotopic, */
                                    /* no pAltInChI->StereoIsotopic or pInChI->StereoIsotopic have been created yet. */
                                    /* added 10-10-2007 to fix i2i bug for Rel/Rac stereo */
                                    if ( bIso && bNo_InChI_t && bNo_InChI_m /* no /t... or /mN in Fixed-H  */ && !nFlags &&
                                         !(pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->t_parity) &&
                                         !(pInChI->StereoIsotopic    && pInChI->StereoIsotopic->t_parity) &&
                                         s[iINChI][jAlt][0]==0 && s[iINChI][j][0]==0 && s[iINChI][j][1] == 0 &&
                                         /* /sN, N=2 or 3 only in Mobile-H isotopic */
                                         (s[iINChI][jAlt][1] & ((INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO))) )
                                    {
                                        /* create two zero/NULL-initialized isotopic stereo if they do not exist */
                                        if ( !pInChI->StereoIsotopic    && 0 > (ret2 = CopySegment( pInChI,    pAltInChI, CPY_SP3_M, 1, -1)) 
                                            /* -- the following will be created later, in TAUT_YES part of the code -- */
                                            /*|| !pAltInChI->StereoIsotopic && 0 > (ret2 = CopySegment( pAltInChI, pAltInChI, CPY_SP3_M, 1, -1))*/ ) {
                                            goto exit_function;
                                        }
                                        /* same value = 2 for MI and FI; here we assign only FI */
                                        pInChI->StereoIsotopic->nCompInv2Abs = 2;
                                        pInChI->Stereo->nCompInv2Abs = 0;
                                        /* -- the following will be assigned later, in TAUT_YES part of the code -- */
                                        /*
                                        pAltInChI->StereoIsotopic->nCompInv2Abs = 2;
                                        pAltInChI->Stereo->nCompInv2Abs = 0;
                                        */
                                    } else
                                        pInChI->Stereo->nCompInv2Abs =  0;
                                    if ( !(pInChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                        pInChI->nFlags |= ((nFlagsAlt|nFlags) & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO));
                                    }
                                }

                            } else
                            /* copy sp3 inversion info: non-isotopic Mobile H => non-isotopic Fixed H  */
                            if ( pAltInChI->Stereo && pAltInChI->Stereo->nNumberOfStereoCenters &&
                                 pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters &&
                                 pAltInChI->Stereo->nCompInv2Abs &&
                                 (!pInChI->Stereo->nCompInv2Abs || NO_VALUE_INT == pInChI->Stereo->nCompInv2Abs) ) {
                                if ( !(nFlagsAlt && !nFlags ) || NO_VALUE_INT == pInChI->Stereo->nCompInv2Abs ) {
                                     /* ??? */
                                    pInChI->Stereo->nCompInv2Abs = pAltInChI->Stereo->nCompInv2Abs;
                                }
                            }

                            /* use same rule to copy stereobonds */
                            if ( pAltInChI->Stereo && pAltInChI->Stereo->nNumberOfStereoBonds &&
                                 (!pInChI->Stereo || !pInChI->Stereo->b_parity) ) {
                                if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP2, 0, 0)) ) {
                                    goto exit_function;
                                }
                            }
                            /* (3) isotopic Mobile H -> isotopic Fixed H */
                            /* if !FH_Stereo && !MH_Stereo && MH_IsoStereo!=NULL && FH_IsoStereo==NULL */
                            if ( bIso ) { 
                                if ( !(pInChI->Stereo && pInChI->Stereo->t_parity) &&                                    /* !FH_Stereo */
                                     !(pAltInChI->Stereo && pAltInChI->Stereo->t_parity) &&                              /* !MH_Stereo */
                                     (pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nNumberOfStereoCenters) && /*  MH_IsoStereo */
                                     (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->t_parity) ) {                  /* !FH_IsoStereo */
                                    /* copy sp3 iso stereo MI->FI (/t) and, if FH nCompInv2Abs (/m) is missing, copy it, too, MI->FI */ 
                                    if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3, 1, 1)) ||
                                         ((!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) &&
                                         0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3_M, 1, 1))) ) {
                                        goto exit_function;
                                    }
                                    /* in case of missing nCompInv2Abs, Relative or Racemic stereo 2005-05-10 */
                                    if ( pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT &&
                                         (nFlagsAlt & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                        pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][jAlt][1]>0? 2 : 0;
                                        if ( !(pInChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                            pInChI->nFlags |= (nFlagsAlt & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO));
                                        }
                                    }
                                } else
                                /* copy sp3 inversion info only: isotopic Mobile H -> isotopic Fixed H  */
                                if ( !(pInChI->Stereo && pInChI->Stereo->t_parity) &&                                     /* !FH_Stereo    /t */
                                     !(pAltInChI->Stereo && pAltInChI->Stereo->t_parity) &&                               /* !MH_Stereo    /t */
                                     (pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nNumberOfStereoCenters) &&  /*  MH_IsoStereo /t */
                                     (pInChI->StereoIsotopic && pInChI->StereoIsotopic->nNumberOfStereoCenters) &&        /*  FH_IsoStereo /t */
                                     pAltInChI->StereoIsotopic->nCompInv2Abs &&                                           /*  MH_IsoStereo /m */
                                     (!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) ) { /*  !FH_IsoStereo /m */
                                    /* added 02-09-2006 */
                                    if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3_M, 1, 1)) ) {
                                        goto exit_function;
                                    }
                                }
                                /* use same rule to copy stereobonds */
                                if ( !(pInChI->Stereo && pInChI->Stereo->b_parity) &&
                                     !(pAltInChI->Stereo && pAltInChI->Stereo->b_parity) &&
                                     (pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nNumberOfStereoBonds) &&
                                     (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->b_parity) ) {
                                    if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP2, 1, 1)) ) {
                                        goto exit_function;
                                    }
                                }

                                /* (4) Copy Fixed-H -> isotopic Fixed-H */
                                /* if FH_Stereo && !MH_IsoStereo && && !FH_IsoStereo */
                                if ( (pInChI->Stereo && pInChI->Stereo->nNumberOfStereoCenters) &&              /* FH_Stereo     /t */
                                     !(pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->t_parity) &&     /* !MH_IsoStereo /t */
                                     !(pInChI->StereoIsotopic && pInChI->StereoIsotopic->t_parity) ) {          /* !FH_IsoStereo /t */

                                    /* added 10-10-2007 DT: copy MH_Iso /m => FH_Iso /m to fix i2i bug for Abs stereo */
                                    /* InChI string contains: MH(/t...), MH_Iso(/mN, no /t), FH(no /t /m), FH_Iso(no /t /m) */
                                    if ( pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nCompInv2Abs && /* MH_IsoStereo /m */
                                         bNo_InChI_t &&
                                         NO_VALUE_INT != pAltInChI->StereoIsotopic->nCompInv2Abs &&              /* undef FH_IsoStereo /m */
                                         !(pInChI->StereoIsotopic && NO_VALUE_INT != pInChI->StereoIsotopic->nCompInv2Abs)) {
                                        if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP3_M, 1, 1))) {
                                            goto exit_function;
                                        }
                                    }

                                    /* added 05-09-2006: copy sp3 FH=>FH_Iso */
                                    if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3, 1, 0)) ||
                                         ((!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) &&
                                         0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3_M, 1, 0))) ) {
                                        goto exit_function;
                                    }
                                    /* in case of missing nCompInv2Abs, Relative or Racemic stereo, /sN in Fixed-H, 2005-05-10 */
                                    if ( pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT &&
                                         (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                        if ( s[iINChI][j][0] > 0 && s[iINChI][j][1] > 0 ) {
                                             /* suppose once in a while only non-taut stereo changes if inverted */
                                            pInChI->StereoIsotopic->nCompInv2Abs = 2;
                                            pInChI->Stereo->nCompInv2Abs =  (++nFH2iFH_AltInv)%2? 2:0;
                                        } else
                                        if ( (s[iINChI][j][0] > 0 || s[iINChI][j][1] > 0) && s[iINChI][j][1] >= 0 ) /* ??? != NO_VALUE_INT ??? */
                                            pInChI->StereoIsotopic->nCompInv2Abs =  2;
                                        else
                                            pInChI->StereoIsotopic->nCompInv2Abs =  0;
                                        if ( !(pInChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) ) {
                                            pInChI->nFlags |= (nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO));
                                        }
                                    }

                                } else
                                /* copy sp3 inversion info only: Fixed-H -> isotopic Fixed H  */
                                if ( (pInChI->Stereo && pInChI->Stereo->t_parity) &&
                                     !(pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->t_parity) &&
                                     (pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nNumberOfStereoCenters) &&
                                     (pInChI->StereoIsotopic && pInChI->StereoIsotopic->nNumberOfStereoCenters) &&
                                     pInChI->Stereo->nCompInv2Abs && 
                                     (!pInChI->StereoIsotopic->nCompInv2Abs || NO_VALUE_INT == pInChI->StereoIsotopic->nCompInv2Abs) ) {
                                    /* added 05-09-2006 */
                                    if ( 0 > (ret2 = CopySegment( pInChI, pInChI, CPY_SP3_M, 1, 0)) ) {
                                        goto exit_function;
                                    }
                                }
                            }
                            if ( bIso && 
                                 !(pInChI->Stereo && pInChI->Stereo->nNumberOfStereoBonds) &&
                                 !(pAltInChI->Stereo && pAltInChI->Stereo->nNumberOfStereoBonds) &&
                                  (pAltInChI->StereoIsotopic && pAltInChI->StereoIsotopic->nNumberOfStereoBonds) &&
                                 (!pInChI->StereoIsotopic || !pInChI->StereoIsotopic->b_parity) ) {
                                if ( 0 > (ret2 = CopySegment( pInChI, pAltInChI, CPY_SP2, 1, 1)) ) {
                                    goto exit_function;
                                }
                            }
                        }
                    }  
                } /* end of component cycle (k) */
            } /* end of Mobile/Fixed H cycle (j) */

            /**** replace NO_VALUE_INT with zeroes in all Mobile & Fixed H components ****/
            for ( j = 0; j < TAUT_NUM; j ++ ) {
                for ( k = 0; k < nNumComponents[iINChI][j]; k ++ ) {
                    if ( pInpInChI[iINChI][j] ) {
                        INChI *pInChI = &pInpInChI[iINChI][j][k];
                        if ( pInChI->nTotalCharge == NO_VALUE_INT ) {
                            pInChI->nTotalCharge = 0;
                        }
                        if ( pInChI->Stereo && pInChI->StereoIsotopic &&
                             pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                                 if ( pInChI->Stereo->nNumberOfStereoCenters &&
                                      pInChI->Stereo->nCompInv2Abs != NO_VALUE_INT ) {
                                 pInChI->StereoIsotopic->nCompInv2Abs = pInChI->Stereo->nCompInv2Abs;
                            }
                        }
                        /* Add special nCompInv2Abs=2 to force /s2 or /s3 in InChI output */
                        if ( pInChI->Stereo && pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT ) {
                            if ( pInChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO) &&
                                 pInChI->Stereo->nNumberOfStereoCenters ) {
                                pInChI->Stereo->nCompInv2Abs = (s[iINChI][j][0]>0 /*|| s[iINChI][j][1]>0*/)? 2 : 0; /* we do not know the real value */
                            } else {
                                pInChI->Stereo->nCompInv2Abs = 0;
                            }
                        }
                        if ( pInChI->StereoIsotopic && pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT ) {
                            if ( pInChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO) &&
                                 pInChI->StereoIsotopic->nNumberOfStereoCenters ) {
                                pInChI->StereoIsotopic->nCompInv2Abs = s[iINChI][j][1]>0? 2 : 0; /* we do not know the real value */
                            } else {
                                pInChI->StereoIsotopic->nCompInv2Abs = 0;
                            }
                        }
                        /* added 02-07-2006 */
                        if ( (pInChI->Stereo && pInChI->Stereo->nCompInv2Abs == NO_VALUE_INT) ||
                             (pInChI->StereoIsotopic && pInChI->StereoIsotopic->nCompInv2Abs == NO_VALUE_INT) ) {
                            ret2 = RI_ERR_PROGR;
                            goto exit_function;
                        }
                        if ( !pInChI->bDeleted && pInChI->nNumberOfAtoms ) {
                            tot_charge[iINChI][j] += pInChI->nTotalCharge;
                            for ( m = 0; m < pInChI->nNumberOfAtoms; m ++ ) {
                                if ( pInChI->nAtom[m] < EL_NUMBER_H || pInChI->nAtom[m] > nElDataLen ) {
                                    ret2 = RI_ERR_PROGR;
                                    goto exit_function;
                                }
                                /* all atoms except H */
                                if ( pInChI->nAtom[m] > EL_NUMBER_H ) {
                                    num_elem[iINChI][j][pInChI->nAtom[m]].num ++;
                                }
                            }
                            if ( 0 > (ret2 = GetInChINumH( pInChI, &m ) ) ) {
                                goto exit_function;
                            }
                            num_elem[iINChI][j][EL_NUMBER_H].num += m;
                        }
                    }
                }
            }


            for ( j = 0; j < TAUT_NUM; j ++ ) {
                for ( k = 0; k < nNumComponents[iINChI][j]; k ++ ) 
                {
                    if ( pInpInChI[iINChI][j] ) 
                    {
                        INChI *pInChI = &pInpInChI[iINChI][j][k];
                        if ( pInChI->Stereo && !pInChI->Stereo->nNumberOfStereoCenters ) 
                        {
                            pInChI->Stereo->nCompInv2Abs = 0;
                        }
                        if ( pInChI->StereoIsotopic && !pInChI->StereoIsotopic->nNumberOfStereoCenters ) 
                        {
                            pInChI->StereoIsotopic->nCompInv2Abs = 0;
                        }
                    }
                }
            }


#if ( FIX_I2I_STEREOCONVERSION_BUG3 == 1 )	
/* (2008-04-10)   1=> Fix bug of i2i conversion */
/* (missed repeating /s in FI after F for multi-component case) */	
            if (nNumComponents[iINChI][TAUT_NON]>1)					 /* if multi-component */
            if ( !s[iINChI][TAUT_YES][0] && !s[iINChI][TAUT_YES][1] )/* if no /s in M, MI */
            if ( (s[iINChI][TAUT_NON][0]>1) && (s[iINChI][TAUT_NON][1]>1) )	 /* if /srel/srac in both F, FI */			
            if ( s[iINChI][TAUT_NON][0] == s[iINChI][TAUT_NON][1] )  /* if same stereo in F and FI */
            /* we assume that at least one component in F has no actual stereo */
            /* and place deliberately 0 to appropriate place */
            for ( k = 0; k < nNumComponents[iINChI][TAUT_NON]; k ++ ) 
            {
                INChI *pInChI = &pInpInChI[iINChI][TAUT_NON][k]; 
                if (pInChI->Stereo->nCompInv2Abs!=0)
                { 
                    pInChI->Stereo->nCompInv2Abs = 0; 
                    goto fini; 
                }
            }
fini:		;
#endif


            if ( num_elem[iINChI][TAUT_YES] ) {
                tot_charge[iINChI][TAUT_YES]  += nNumProtons[iINChI][TAUT_YES].nNumRemovedProtons;
                num_elem[iINChI][TAUT_YES][EL_NUMBER_H].num += nNumProtons[iINChI][TAUT_YES].nNumRemovedProtons;
            }

            /**** Count H and isotopic H in Mobile and Fixed H represntations of components */
            /* if at least one component has Fixed-H layer then all components have Fixed-H */
            /* layer; those whose Fixed-H layer is empty have Fixed-H layer same as Mobile-H layer */
            if ( nNumComponents[iINChI][TAUT_NON] ) {
                /* only if both Mobile and Fixed H exist */
                int nFormulaH[TAUT_NUM], nNumH[TAUT_NUM], nCharge[TAUT_NUM], nNumIsotopicH[TAUT_NUM][NUM_H_ISOTOPES];
                int nRemovedCharge, nRemovedH, nRemovedIsotopicH[NUM_H_ISOTOPES], nFoundRemovedIsoH;
                int nTotRemovedProtons, nTotRemovedIsotopicH[NUM_H_ISOTOPES], bExists[TAUT_NUM];
                INChI *pInChI[TAUT_NUM];
                nTotRemovedProtons = 0;
                memset( nTotRemovedIsotopicH, 0, sizeof(nTotRemovedIsotopicH) );
                len2 = inchi_max( nNumComponents[iINChI][TAUT_YES], nNumComponents[iINChI][TAUT_NON] );
                
                for ( k = 0; k < len2; k ++ ) {
                    /* k is a component index */
                    for ( j = 0; j < TAUT_NUM; j ++ ) {
                        /* j is 0=TAUT_NON or 1=TAUT_YES */
                        pInChI[j]  = NULL; /* initialization 2006-03 */
                        bExists[j] = (k < nNumComponents[iINChI][j]) &&
                                     pInpInChI[iINChI][j][k].nNumberOfAtoms &&
                                    !pInpInChI[iINChI][j][k].bDeleted;
                    }
                    if ( !bExists[TAUT_NON] ) {
                        /* TAUT_YES does not exist for a proton (H+) in TAUT_NON */
                        ret2 = RI_ERR_SYNTAX;
                        goto exit_function;
                    }
                    /* at this point at least one of Mobile[k] and Fixed[k] real InChI exists */
                    /* initialize for counting removed protons and isotopic H from kth  Mobile-H component */
                    for ( j = 0; j < TAUT_NUM; j ++ ) {
                        if ( bExists[j] ) {
                            pInChI[j] = &pInpInChI[iINChI][j][k];  /* BC: reading uninit memory (fixed?) */
                        }
                        nFormulaH[j] = 0;
                        nNumH[j]     = 0;
                        nCharge[j]   = 0;
                        for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                            nNumIsotopicH[j][m] = 0;
                        }
                    }
                    /* extract number of H, isotopic H, and charge */
                    for ( j = 0; j < TAUT_NUM; j ++ ) {
                        if ( !bExists[j] )
                            continue;
                        if ( 0 > (ret2 = GetInChIFormulaNumH( pInChI[j], &nFormulaH[j] )) ||
                             0 > (ret2 = GetInChINumH( pInChI[j], &nNumH[j] )) ||
                             0 > (ret2 = GetInChIIsoH( pInChI[j], nNumIsotopicH[j] )) ) {
                            goto exit_function;
                        }
                        nCharge[j] = pInChI[j]->nTotalCharge;
                    }
                    for ( j = 0; j < TAUT_NUM; j ++ ) {
                        if ( !bExists[j] )
                            continue;
                        if ( nFormulaH[j] != nNumH[j] ) {
                            ret2 = RI_ERR_SYNTAX;
                            goto exit_function;
                         }
                    }
                    nFoundRemovedIsoH = 0;
                    nRemovedCharge = nCharge[TAUT_NON] - nCharge[TAUT_YES];
                    nRemovedH      = nNumH[TAUT_NON]   - nNumH[TAUT_YES];
                    for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                        nFoundRemovedIsoH += 0 != (nRemovedIsotopicH[m] = nNumIsotopicH[TAUT_NON][m] -
                                                                          nNumIsotopicH[TAUT_YES][m] );
                    }
                    if ( nRemovedCharge != nRemovedH ) {
                        ret2 = RI_ERR_SYNTAX;
                        goto exit_function;
                    }
                    if ( nRemovedCharge || nFoundRemovedIsoH ) {
                        COMPONENT_REM_PROTONS *pNumProtons;
                        if ( !nNumProtons[iINChI][TAUT_YES].pNumProtons ) {
                            /* allocate only if needed */
                            nNumProtons[iINChI][TAUT_YES].pNumProtons = 
                                (COMPONENT_REM_PROTONS *) inchi_calloc(len2,
                                                     sizeof(nNumProtons[0][0].pNumProtons[0]));
                            if ( !nNumProtons[iINChI][TAUT_YES].pNumProtons ) {
                                ret2 = RI_ERR_ALLOC;
                                goto exit_function;
                            }
                        }
                        pNumProtons = nNumProtons[iINChI][TAUT_YES].pNumProtons+k;
                        pNumProtons->nNumRemovedProtons = nRemovedH;
                        nTotRemovedProtons += nRemovedH;
                        for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                            pNumProtons->nNumRemovedIsotopicH[m] = nRemovedIsotopicH[m];
                            nTotRemovedIsotopicH[m] += nRemovedIsotopicH[m];
                        }
                        /* make sure the Mobile-H InChI has nTautomer */
                        if ( pInChI[TAUT_YES] && bExists[TAUT_YES] ) {
                            if ( !pInChI[TAUT_YES]->lenTautomer ) {
                                pInChI[TAUT_YES]->lenTautomer = 1;
                            }
                            if ( !pInChI[TAUT_YES]->nTautomer ) {
                                pInChI[TAUT_YES]->nTautomer = (AT_NUMB *)inchi_calloc(pInChI[TAUT_YES]->lenTautomer, sizeof(pInChI[0]->nTautomer[0]) );
                            }
                        }
                    }
                }
                if ( nNumProtons[iINChI][TAUT_YES].pNumProtons ) {
                    /* check consistency */
#if ( FIX_ISO_FIXEDH_BUG_READ == 1 )
                    int iso_diff[NUM_H_ISOTOPES], iso_diff_tot=0;
#endif
                    if ( nTotRemovedProtons != nNumProtons[iINChI][TAUT_YES].nNumRemovedProtons ) {
                        ret2 = RI_ERR_SYNTAX;
                        goto exit_function;
                    }
#if ( FIX_ISO_FIXEDH_BUG_READ == 1 )
                    for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                        iso_diff[m] = nNumProtons[iINChI][TAUT_YES].nNumRemovedIsotopicH[m]-nTotRemovedIsotopicH[m];
                        if ( iso_diff[m] < 0 )
                        {
                            ret2 = RI_ERR_SYNTAX;
                            goto exit_function;
                        } else {
                            /* InChI-1.02b bug: nTotRemovedIsotopicH[m] < nNumProtons[iINChI][TAUT_YES].nNumRemovedIsotopicH[m] */
                            /* in non-tautomeric components where D(+) or T(+) was removed from -NH(+)= or =OH(+)           */ 
                            iso_diff_tot += iso_diff[m]; 
                        }
                    }
                    if ( iso_diff_tot ) {
                        if ( 0 > bIsoMayBeArranged( bInchi2Struct, iso_diff, nNumProtons, pInpInChI, nNumComponents, iINChI )) {
                            ret2 = RI_ERR_SYNTAX;
                            goto exit_function;
                        }
                    }
#else
                    for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                        if ( nTotRemovedIsotopicH[m] != nNumProtons[iINChI][TAUT_YES].nNumRemovedIsotopicH[m] ) {
                            ret2 = RI_ERR_SYNTAX;
                            goto exit_function;
                        }
                    }
#endif
                }
            }

            /* make Mobile H and Fixed H InChI arrays have same length */
            len2 = len1 = 0;
            if ( nNumComponents[iINChI][TAUT_YES] < nNumComponents[iINChI][TAUT_NON] ) {
                j = TAUT_YES; /* less components in Mobile-H layer */
                len2 = nNumComponents[iINChI][TAUT_NON];
                len1 = nNumComponents[iINChI][TAUT_YES];
            } else
            if ( nNumComponents[iINChI][TAUT_YES] > nNumComponents[iINChI][TAUT_NON] ) {
                j = TAUT_NON;  /* less components in Fixed-H layer */
                len2 = nNumComponents[iINChI][TAUT_YES];
                len1 = nNumComponents[iINChI][TAUT_NON];
            }
            /* always len1 <= len2; if Mobile-H and Fixed-H have same number of components then len1=len2=0  */
            if ( len2 && len1 ) {
                INChI *pInChI = (INChI *)inchi_calloc( len2, sizeof( pInChI[0] ) );
                if ( !pInChI ) {
                    ret2 = RI_ERR_ALLOC;
                    goto exit_function;
                }
                memcpy( pInChI, pInpInChI[iINChI][j], len1 * sizeof( pInChI[0] ) );
                inchi_free( pInpInChI[iINChI][j] );
                pInpInChI[iINChI][j] = pInChI;
                nNumComponents[iINChI][j] = len2;
                for ( ; len1 < len2; len1 ++ ) {
                    if ( j == TAUT_YES ) {
                        /* mark added to Mobile H layer components as deleted protons */
                        if ( 0 > (ret2 = nFillOutProtonMobileH( pInpInChI[iINChI][j]+len1 ) ) ) {
                            goto exit_function;
                        }
                        if ( 0 > (ret2 = nProtonCopyIsotopicInfo( pInpInChI[iINChI][j]+len1/* to */,
                                                                  pInpInChI[iINChI][TAUT_NON]+len1/* from */ ) ) ) {
                            goto exit_function;
                        }
                    } else {
                        /* mark added to Fixed H layer components as empty deleted */
                        /* this should not happen */
                        pInChI[len1].bDeleted = 1;
                    }
                }
            }
        } /* end of iINChI cycle */
        /* check balances */
        for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) {
            for ( i = iINChI; i < INCHI_NUM; i ++ ) {
                for ( j = 0; j < TAUT_NUM; j ++ ) {
                    for ( k = j; k < TAUT_NUM; k ++ ) {
                        if ( (iINChI!=i || j !=k) && num_elem[iINChI][j] && num_elem[i][k] ) {
                            if ( tot_charge[iINChI][j] != tot_charge[i][k] ) {
                                ret2 = RI_ERR_SYNTAX;
                                goto exit_function;
                            }
                            for ( m = 0; m <= nElDataLen; m ++ ) {
                                if ( num_elem[iINChI][j][m].num != num_elem[i][k][m].num ) {
                                    ret2 = RI_ERR_SYNTAX;
                                    goto exit_function;
                                }
                            }
                            /*
                            if ( memcmp( num_elem[iINChI], num_elem[i][k], (nElDataLen+1)*sizeof(num_elem[0][0][0]) ) {
                                ret2 = RI_ERR_SYNTAX;
                                goto exit_function;
                            }
                            */
                        }
                    }
                }
            }
        }


    } else {
        ret2 = ret;
    }
exit_function:
    for ( i = 0; i < INCHI_NUM; i ++ ) {
        for ( j = 0; j < TAUT_NUM; j ++ ) {
            if ( num_elem[i][j] ) {
                inchi_free( num_elem[i][j] );
                num_elem[i][j] = NULL;
            }
        }
    }
    *nErr = (ret2 < 0 && ret2 != RI_ERR_EOL)? ret2 : 0;
    return ret;
}

/**************************************************************************************/
#if ( FIX_ISO_FIXEDH_BUG_READ == 1 )
#undef TAUT_YES
int bIsoMayBeArranged( int bInchi2Struct, int iso_diff[NUM_H_ISOTOPES], REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                       INChI *pInpInChI[INCHI_NUM][TAUT_NUM], int nNumComponents[INCHI_NUM][TAUT_NUM], int iINChI )
{
const int TAUT_YES = 1;
    int i, k, m, n_found=0, n_found_at_in_component, n_found_H_in_component, i_iso_at, num_iso_H=0, num_iso_H_orig, num_add_iso_H, orig_add_H;
    for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
        num_iso_H += iso_diff[m];
    }
    num_iso_H_orig = num_iso_H;
    for ( k = 0; k < nNumComponents[iINChI][TAUT_YES] && k < nNumComponents[iINChI][TAUT_NON]; k ++ ) {
        INChI *pInChI = &pInpInChI[iINChI][TAUT_NON][k];
        INChI *pInChITaut = &pInpInChI[iINChI][TAUT_YES][k];
        if ( pInChITaut->bDeleted || pInChI->bDeleted ||
             pInChITaut->nNumberOfIsotopicAtoms > 0 ||
             pInChITaut->lenTautomer > 1 && pInChITaut->nTautomer && pInChITaut->nTautomer[0] > 0 ||
             NULL == nNumProtons[iINChI][TAUT_YES].pNumProtons ||
             nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons <= 0 ||
             pInChI->nNumberOfIsotopicAtoms > 0 ||
             nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[0] ||
             nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[1] ||
             nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[2]
             ) {
            continue;
        }
        /* check if fixed-H has isotopic H; count the possibilities */
        orig_add_H = nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons;
        n_found_at_in_component = 0; /* number of atoms that may accept isotopic H */
        n_found_H_in_component = 0;
        for ( i = 0; i < pInChI->nNumberOfAtoms; i ++ ) {
            int nNumRemovedH = (int)pInChI->nNum_H[i] - (int)pInChITaut->nNum_H[i];
            if ( nNumRemovedH > 0) {
                n_found_at_in_component ++;
                n_found_H_in_component += nNumRemovedH;
            }
        }
        if ( n_found_at_in_component > 0 && num_iso_H > 0 && bInchi2Struct ) {
            pInChI->IsotopicAtom = (INChI_IsotopicAtom *)calloc(inchi_min(n_found_at_in_component, num_iso_H), sizeof(pInChI->IsotopicAtom[0]));
        }
        for ( i = 0, i_iso_at = 0; i < pInChI->nNumberOfAtoms; i ++ ) {
            int nNumRemovedH = (int)pInChI->nNum_H[i] - (int)pInChITaut->nNum_H[i];
            n_found += nNumRemovedH; /* found H removed in mobile-H layer */
            if ( nNumRemovedH > 0 && num_iso_H > 0 && orig_add_H ) {
                for ( m = 0; m < NUM_H_ISOTOPES && 0 < num_iso_H && 0 < orig_add_H && 0 < nNumRemovedH; m ++ ) {
                    if ( iso_diff[m] > 0 ) {
                        num_add_iso_H = inchi_min( iso_diff[m], nNumRemovedH ); /* atom limit */
                        if ( num_add_iso_H > orig_add_H )                       /* component limit */
                            num_add_iso_H = orig_add_H;
                        iso_diff[m] -= num_add_iso_H;                           /* update tot removed single isotope H limit */
                        num_iso_H   -= num_add_iso_H;                           /* update tot removed isotopic H limit */
                        orig_add_H  -= num_add_iso_H;                           /* update component limit */
                        nNumRemovedH -= num_add_iso_H;                          /* update atom limit */
                        nNumProtons[iINChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[m] += num_add_iso_H;
                        if ( pInChI->IsotopicAtom ) {
                            pInChI->IsotopicAtom[i_iso_at].nAtomNumber = i+1;
                            switch( m ) {
                            case 0:
                                pInChI->IsotopicAtom[i_iso_at].nNum_H += num_add_iso_H;
                                break;
                            case 1:
                                pInChI->IsotopicAtom[i_iso_at].nNum_D += num_add_iso_H;
                                break;
                            case 2:
                                pInChI->IsotopicAtom[i_iso_at].nNum_T += num_add_iso_H;
                                break;
                            }
                        }
                    }
                }
                if ( pInChI->IsotopicAtom ) {
                    i_iso_at ++;
                }
            }
        }
        if ( pInChI->IsotopicAtom && i_iso_at ) {
            pInChI->nNumberOfIsotopicAtoms = i_iso_at;
        }
    }
    if ( n_found - num_iso_H >= 0 ) {
        /* Success. Arrange isotopic H between components */

    }

    return n_found - num_iso_H_orig; /* >0 => ambiguous reconstruction, 0 => unambiguous, <0 => impossible */
}
#define TAUT_YES 1
#endif

/******************************************************************************************************/
typedef enum tagAuxInfoState {
    AST_VERSION,                         /* 0    */
    
    AST_MOBILE_H_NUMBERS,                /* 1 /N:  */
    AST_MOBILE_H_ATOM_EQ,                /* 2 /E:  */
    AST_MOBILE_H_GROUP_EQ,               /* 3 /gE: */
    AST_MOBILE_H_SP3_INV,                /* 4 /it: */
    AST_MOBILE_H_SP3_INV_NUMBERS,        /* 5 /iN: */    
    
    AST_MOBILE_H_ISO_LAYER_FORK,         /* 6 */

    AST_MOBILE_H_ISO_NUMBERS,            /* 7 /I:  */
    AST_MOBILE_H_ISO_ATOM_EQ,            /* 8 /E:  */
    AST_MOBILE_H_ISO_GROUP_EQ,           /* 9 /gE: */
    AST_MOBILE_H_ISO_SP3_INV,            /* 10 /it: */
    AST_MOBILE_H_ISO_SP3_INV_NUMBERS,    /* 11 /iN: */

    AST_FIXED_H_LAYER_FORK,              /* 12 */

    AST_FIXED_H_NUMBERS,                 /* 13 /F:  */
    AST_FIXED_H_ATOM_EQ,                 /* 14 /E:  */
    AST_FIXED_H_SP3_INV,                 /* 15 /it: */
    AST_FIXED_H_SP3_INV_NUMBERS,         /* 16 /iN: */

    AST_FIXED_H_ISO_LAYER_FORK,          /* 17 */

    AST_FIXED_H_ISO_NUMBERS,             /* 18 /I:  */
    AST_FIXED_H_ISO_ATOM_EQ,             /* 19 /E:  */
    AST_FIXED_H_ISO_SP3_INV,             /* 20 /it: */
    AST_FIXED_H_ISO_SP3_INV_NUMBERS,     /* 21 /iN: */

    AST_REVERSE_INFO_CRV,                      /* 22 /CRV: */
    AST_REVERSE_INFO_ATOMS,                    /* 23 /rA:  */
    AST_REVERSE_INFO_BONDS,                    /* 24 /rB:  */
    AST_REVERSE_INFO_XYZ,                      /* 25 /rC:  */

    AST_RECONNECTED_LAYER_FORK,          /* 26 /R:   */
    AST_RECONNECTED_LAYER_NUMBERS        /* 27 */

}AUX_INFO_STATE;

/************************************************************************************/
int ParseAuxSegmentVersion( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    const char *q;
    if ( isdigit( UCINT *str ) && (inchi_strtol( str, &q, 10), !*q) ) {
        return 1;
    }
    return RI_ERR_SYNTAX;
}
/************************************************************************************/
int CopyAtomNumbers( INChI *pInChI_To, int bIsoTo, INChI *pInChI_From, int bIsoFrom )
{
    AT_NUMB *pTo, *pFrom;
    if ( !pInChI_To || !pInChI_From || pInChI_To->bDeleted || pInChI_From->bDeleted ||
         !pInChI_To->nNumberOfAtoms || !pInChI_From->nNumberOfAtoms ||
         pInChI_To->nNumberOfAtoms  != pInChI_From->nNumberOfAtoms ||
         !pInChI_From->nPossibleLocationsOfIsotopicH ) {
        return RI_ERR_PROGR;
    }
    if ( !pInChI_To->nPossibleLocationsOfIsotopicH ) {
        pInChI_To->nPossibleLocationsOfIsotopicH = (AT_NUMB *)inchi_calloc( 2*pInChI_To->nNumberOfAtoms,
                                                   sizeof(pInChI_To->nPossibleLocationsOfIsotopicH[0]));
        if ( !pInChI_To->nPossibleLocationsOfIsotopicH ) {
            return RI_ERR_ALLOC;
        }
    }
    pTo   = pInChI_To->nPossibleLocationsOfIsotopicH + (bIsoTo? 0 : pInChI_To->nNumberOfAtoms );
    pFrom = pInChI_From->nPossibleLocationsOfIsotopicH + (bIsoFrom? 0 : pInChI_To->nNumberOfAtoms );
    if ( pTo == pFrom ) {
        return RI_ERR_PROGR;
    }
    memcpy( pTo, pFrom, pInChI_To->nNumberOfAtoms*sizeof(pTo[0]) );
    return 1;
}
/************************************************************************************/
int ParseAuxSegmentNumbers( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state, int *pbAbc )
{
    int bIso = 0, iComponent = 0, nNumComponents, bIso_From, bAltInChIExists;
    INChI   *pInChI = NULL, *pAltInChI = NULL, *pInChI_From = NULL;
    const char  *p, *q, *pStart, *pEnd, *t;
    static const char  mult_type[] = "mnM";
    int      val, ret, k, mpy_component, num;
    AT_NUMB *pNumb;
    int      base = 10;
    /* save isotopic numbering into the first nNumberOfAtoms elements of INChI::nPossibleLocationsOfIsotopicH */
    /* save non-isotopic numbering into the second half of nNumberOfAtoms elements of INChI::nPossibleLocationsOfIsotopicH */

    switch( state ) {
    case AST_MOBILE_H_NUMBERS:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "N:", 2 ) )
            return 0;
        break;
    case AST_FIXED_H_NUMBERS:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "F:", 2 ) )
            return 0;
        break;
    case AST_MOBILE_H_ISO_NUMBERS:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "I:", 2 ) )
            return 0;
        bIso = 1;
        break;
    case AST_FIXED_H_ISO_NUMBERS:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "I:", 2 ) )
            return 0;
        bIso = 1;
        break;
    default:
        return RI_ERR_PROGR;
    }
    pStart = str+2;
    if ( !*pStart ) {
        return 1;
    }
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    bAltInChIExists = (NULL != pInpInChI[ALT_TAUT(bMobileH)]);
    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        /* check */
        if ( !pInpInChI[bMobileH] ) {
            return 1; /* invalid aux info */
        }
        pInChI    = pInpInChI[bMobileH] + iComponent;
        pAltInChI = pInpInChI[ALT_TAUT(bMobileH)] + iComponent;
        if ( ((isdigit(UCINT *pStart) &&
             0 < (val = (int)inchi_strtol( pStart, &q, 10))) ||
             (q = pStart, val=1) )&&
             (t=strchr(mult_type, *q)) && q+1 == pEnd ) {
            
            /* process the abbreviation */
            
            pInChI_From = NULL;

            switch( bMobileH ) {

            case TAUT_YES:
                switch ( bIso ) {
                case 0:
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                case 1:
                    if ( *q != 'm' ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    /* isotopic Mobile-H  <-- non-isotopic Mobile H  */
                    pInChI_From = pInChI;
                    bIso_From   = 0;
                    break;
                default:
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
                break;

            case TAUT_NON:
                switch ( *q ) {
                case 'm':  /* same as mobile H */
                    switch( bIso ) {
                    case 0: /* from Mobile-H not isotopic */
                        pInChI_From = bAltInChIExists? pAltInChI:NULL;
                        bIso_From   = 0;
                        break;

                    case 1:
                        pInChI_From = bAltInChIExists? pAltInChI:NULL;;
                        bIso_From   = 1;
                        break;
                    default:
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    break;
                case 'n': /* same as non-isotopic Fixed-H */
                    switch ( bIso ) {
                    case 0:
                        ret = 1; /*RI_ERR_SYNTAX;*/
                        goto exit_function;
                    case 1:
                        pInChI_From = pInChI;
                        bIso_From   = 0;
                    default:
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    break;
                case 'M':  /* same as isotopic Mobile-H */
                    switch ( bIso ) {
                    case 0:
                        ret = RI_ERR_SYNTAX;
                        goto exit_function;
                    case 1:
                        pInChI_From = bAltInChIExists? pAltInChI:NULL;;
                        bIso_From   = 1;
                        break;
                    default:
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    break;
                default:
                    ret = 1; /*RI_ERR_SYNTAX;*/
                    goto exit_function;
                }
                break;

            }
            /* copy */
            if ( pInChI_From ) {
                for ( k = 0; k < val; k ++ ) {
                    CopyAtomNumbers( pInChI+k, bIso, pInChI_From+k, bIso_From );
                }
            }
            mpy_component = val;
        } else {
            mpy_component = 1;
            p = pStart;
            pNumb = pInChI->nPossibleLocationsOfIsotopicH;
            if ( !pNumb ) {
                pNumb = (AT_NUMB *)inchi_calloc( 2*pInChI->nNumberOfAtoms, sizeof(pNumb[0] ) );
                if ( !pNumb ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                pInChI->nPossibleLocationsOfIsotopicH = pNumb;
            }
            pNumb += bIso? 0 : pInChI->nNumberOfAtoms;
            if ( pStart < pEnd && *pbAbc == -1 ) {
                /* check if compressed InChI */
                *pbAbc = isupper( UCINT *pStart)? 1 : 0;
            }
            base = (*pbAbc==1)? ALPHA_BASE : 10;

            if ( *pbAbc == 1 ) {
                for ( k = 0, p = pStart; k < pInChI->nNumberOfAtoms && p < pEnd; k ++, p ++ ) {
                    num = (AT_NUMB)inchi_strtol( p, &q, base );
                    if ( num <= 0 || p == q ) {
                        ret = RI_ERR_SYNTAX;
                        goto exit_function;
                    }
                    pNumb[k] = (AT_NUMB)num;
                    p = q;
                    if ( p == pEnd ) {
                        break; /* main end of cycle */
                    }
                }
            } else {
                for ( k = 0, p = pStart; k < pInChI->nNumberOfAtoms && p < pEnd; k ++, p ++ ) {
                    pNumb[k] = (AT_NUMB)inchi_strtol( p, &q, 10 );
                    p = q;
                    if ( p == pEnd ) {
                        break; /* main end of cycle */
                    } else
                    if ( *p != ',' ) {
                        ret = RI_ERR_SYNTAX;
                        goto exit_function;
                    }
                }
            }
            if ( p != pEnd || k+1 != pInChI->nNumberOfAtoms ) {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
            }
        }

        iComponent += mpy_component;
        if ( *pEnd ) {
            pStart = pEnd+1;
            continue;
        } else {
            break;
        }
    }

    if ( nNumComponents != iComponent ) {
        ret = 1; /*RI_ERR_SYNTAX;*/ /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    return ret;   

}
/**********************************************************************************************************/
int ParseAuxSegmentAtomEqu( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_MOBILE_H_ATOM_EQ:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "E:", 2 ) )
            return 0;
        break;
    case AST_MOBILE_H_ISO_ATOM_EQ:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "E:", 2 ) )
            return 0;
        break;
    case AST_FIXED_H_ATOM_EQ:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "E:", 2 ) )
            return 0;
        break;
    case AST_FIXED_H_ISO_ATOM_EQ:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "E:", 2 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentGroupEqu( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_MOBILE_H_GROUP_EQ:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "gE:", 3 ) )
            return 0;
        break;
    case AST_MOBILE_H_ISO_GROUP_EQ:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "gE:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentSp3Inv( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_MOBILE_H_SP3_INV:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "it:", 3 ) )
            return 0;
        break;
    case AST_MOBILE_H_ISO_SP3_INV:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "it:", 3 ) )
            return 0;
        break;
    case AST_FIXED_H_SP3_INV:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "it:", 3 ) )
            return 0;
        break;
    case AST_FIXED_H_ISO_SP3_INV:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "it:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentSp3InvNumbers( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_MOBILE_H_SP3_INV_NUMBERS:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "iN:", 3 ) )
            return 0;
        break;
    case AST_MOBILE_H_ISO_SP3_INV_NUMBERS:
        if ( bMobileH != TAUT_YES )
            return RI_ERR_PROGR;
        if ( memcmp( str, "iN:", 3 ) )
            return 0;
        break;
    case AST_FIXED_H_SP3_INV_NUMBERS:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "iN:", 3 ) )
            return 0;
        break;
    case AST_FIXED_H_ISO_SP3_INV_NUMBERS:
        if ( bMobileH != TAUT_NON )
            return RI_ERR_PROGR;
        if ( memcmp( str, "iN:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentReverseCRV( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_REVERSE_INFO_CRV:
        if ( memcmp( str, "CRV:", 4 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentReverseAtoms( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_REVERSE_INFO_ATOMS:
        if ( memcmp( str, "rA:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentReverseBonds( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    switch( state ) {
    case AST_REVERSE_INFO_BONDS:
        if ( memcmp( str, "rB:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    return 1;
}
/***********************************************************************************************************/
int ParseAuxSegmentReverseXYZ( const char *str, int bMobileH, XYZ_COORD **ppXYZ, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    const char *pStart, *p, *q;
    XYZ_COORD *pXYZ = NULL;
    int     nLenXYZ=0, i, j;
    switch( state ) {
    case AST_REVERSE_INFO_XYZ:
        if ( memcmp( str, "rC:", 3 ) )
            return 0;
        break;
    default:
        return RI_ERR_PROGR;
    }
    pStart = str+3;
    /* count coordinates */
    for ( p = pStart, nLenXYZ = 0; *p; p ++ ) {
        nLenXYZ += ( *p == ';' );
    }
    if ( !nLenXYZ ) {
        return RI_ERR_SYNTAX;
    }
    if ( NULL == (pXYZ = (XYZ_COORD *)inchi_calloc( nLenXYZ, sizeof(pXYZ[0]) )) ) {
        return RI_ERR_ALLOC;
    }
    for ( p = pStart, i = 0; *p && i < nLenXYZ; p ++, i ++ ) {
        for ( j = 0; j < 3; j ++ ) {
            pXYZ[i].xyz[j] = inchi_strtod( p, &q );
            p = q + (*q == ',' );
        }
        if ( *p != ';' ) {
            break;
        }
    }
    if ( i != nLenXYZ || *p ) {
        return RI_ERR_SYNTAX;
    }
    *ppXYZ = pXYZ;
    return nLenXYZ+1;
}
/************************************************************************************/
int AddAuxSegmentCoord( int nRet, XYZ_COORD *pXYZ, int nLenXYZ, INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                        int nNumComponents[INCHI_NUM][TAUT_NUM] )
{
    int iINChI, j, k, n, m, numAt[TAUT_NUM], num_at, nNumMissingNumbers = 0, ret = 0;
    INChI *pInChI    = NULL;
    INChI *pAltInChI = NULL;
    XYZ_COORD *pxyz;

    /* propagate numberings */
    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) {
        for ( j = TAUT_YES; TAUT_NON <= j; j -- ) {
            for ( k = 0; k < nNumComponents[iINChI][j]; k ++ ) {
                int   jj = ALT_TAUT(j);
                pInChI    = pInpInChI[iINChI][j] + k;
                pAltInChI = (k < nNumComponents[iINChI][jj])? pInpInChI[iINChI][jj] + k : NULL;
                numAt[j]  = ( !pInChI->bDeleted )? pInChI->nNumberOfAtoms : 0;
                numAt[jj] = ( pAltInChI && !pAltInChI->bDeleted )? pAltInChI->nNumberOfAtoms : 0;
                switch( j ) {
                case TAUT_YES:
                    if ( !numAt[j] ) {
                        break; /* component does not exist */
                    }
                    if ( !pInChI->nPossibleLocationsOfIsotopicH ) {
                        nNumMissingNumbers ++;
                        break;
                    }
                    if ( !pInChI->nPossibleLocationsOfIsotopicH[0] ) {
                        if ( pInChI->nPossibleLocationsOfIsotopicH[numAt[j]] ) {
                            /* copy from non-isotopic (2nd half of the at. numbers array) to the isotopic (1st half) */
                            ret = CopyAtomNumbers( pInChI, 1, pInChI, 0 );
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                        } else {
                            inchi_free( pInChI->nPossibleLocationsOfIsotopicH );
                            pInChI->nPossibleLocationsOfIsotopicH = NULL;
                            nNumMissingNumbers ++;
                        }
                    }
                    break;

                case TAUT_NON:
                    if ( !numAt[j] ) {
                        break; /* component does not exist */
                    }
                    if ( !pInChI->nPossibleLocationsOfIsotopicH ) {
                        /* trying to get numbers from Mobile-H component */
                        if ( !numAt[jj] || !(pAltInChI->nPossibleLocationsOfIsotopicH) ) {
                            nNumMissingNumbers ++;
                            break;
                        }
                        if ( pAltInChI->nPossibleLocationsOfIsotopicH[0] ) {
                            ret = CopyAtomNumbers( pInChI, 1, pAltInChI, 1 );
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                        } else
                        if ( pAltInChI->nPossibleLocationsOfIsotopicH[numAt[jj]] ) {
                            ret = CopyAtomNumbers( pInChI, 1, pAltInChI, 0 );
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                        } else {
                            /* pAltInChI->nPossibleLocationsOfIsotopicH should have */
                            /* been deallocated on previous TAUT_YES pass           */
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    } else
                    if ( !pInChI->nPossibleLocationsOfIsotopicH[0] ) {
                        if ( pInChI->nPossibleLocationsOfIsotopicH[numAt[j]] ) {
                            /* copy from non-isotopic to isotopic */
                            ret = CopyAtomNumbers( pInChI, 1, pInChI, 0 );
                            if ( ret < 0 ) {
                                goto exit_function;
                            }
                        } else {
                            inchi_free( pInChI->nPossibleLocationsOfIsotopicH );
                            pInChI->nPossibleLocationsOfIsotopicH = NULL;
                            nNumMissingNumbers ++;
                        }
                    }
                    break;
                }
            }
        }
    }
    /* add coordinates */
    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) {
        for ( j = 0; j < TAUT_NUM; j ++ ) {
            for ( k = 0; k < nNumComponents[iINChI][j]; k ++ ) {
                pInChI    = pInpInChI[iINChI][j] + k;
                num_at    = ( !pInChI->bDeleted )? pInChI->nNumberOfAtoms : 0;
                if ( !num_at ) {
                    if ( pInChI->nPossibleLocationsOfIsotopicH ) {
                        inchi_free( pInChI->nPossibleLocationsOfIsotopicH );
                        pInChI->nPossibleLocationsOfIsotopicH = NULL;
                    }
                    continue;
                }
                if ( !pInChI->nPossibleLocationsOfIsotopicH ) {
                    continue;
                }
                if ( iINChI == INCHI_BAS && num_at == 1 &&
                     pInChI->szHillFormula && !strcmp(pInChI->szHillFormula, "H") &&
                     (int)pInChI->nPossibleLocationsOfIsotopicH[0]-1 >= nLenXYZ ) {
                    ; /* a single atom H disconnected from a metal atom has no coordinates */
                } else {
                    /* add atom coordinates */
                    pxyz = (XYZ_COORD *)inchi_calloc( num_at, sizeof(pxyz[0]));
                    if ( !pxyz ) {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                    for ( n = 0; n < num_at; n ++ ) {
                        m = (int)pInChI->nPossibleLocationsOfIsotopicH[n]-1;
                        if ( m < 0 || m >= nLenXYZ ) {
                            inchi_free( pxyz );
                            ret = RI_ERR_SYNTAX;
                            goto exit_function;
                        }
                        pxyz[n] = pXYZ[m];
                    }
                    pInChI->IsotopicTGroup = (INChI_IsotopicTGroup *)pxyz;
                }
                inchi_free( pInChI->nPossibleLocationsOfIsotopicH );
                pInChI->nPossibleLocationsOfIsotopicH = NULL;
            }
        }
    }
    ret = nRet; /* normal exit */

exit_function:
    return ret;
}
/************************************************************************************/
int ReadInChICoord( INCHI_IOSTREAM *pInp, 
                    SEGM_LINE *pLine, int *pState,
                    INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                    int nNumComponents[INCHI_NUM][TAUT_NUM] )
{
    int   c, fst, ret=RI_ERR_ALLOC;
    int   bMobileH = TAUT_YES, bReconn = INCHI_BAS;
    const char szToken[] = INCHI_TOKEN;
    int  state=-1, prev_state=-1;
    XYZ_COORD *pXYZ = NULL;
    int       nLenXYZ = 0;
    int       bAbc    = -1; /* initially undefined */

    *pState = 0;

    INCHI_HEAPCHK
    /* Get "InChI=1/" */
    if ( pLine->len ) {
        c = pLine->c;
    } else {
        c = nGetInChISegment( pInp, pLine, szToken );
    }
    if ( c == RI_ERR_EOF && !pLine->len && !pLine->str[0] ) {
        ret = c;
        pLine->len = 0;
        goto exit_error;
    }
    if ( pLine->len == 0 || (c != SEG_END && c != RI_ERR_EOF && !INCHI_INP_EOL(c)) ) {
        *pState = -1;
        pLine->len = 0;
        ret = RI_ERR_PROGR;
        goto exit_error;
    }
    if ( memcmp(pLine->str, "AuxInfo=", 8) ) {
        *pState = -1;
        return c;
    }
    state = AST_VERSION;
    ret  = 1; /* means read the next segment */
    do {
        /* read the next segment up to the '/' */
        INCHI_HEAPCHK
        if ( ret < 0 ) {
            *pState = prev_state;
            break;
        }
        prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
        /* prev_state = state;*/
        if ( 0 < ret ) {
            /* read next segment */
            if ( c != RI_ERR_EOF && c != SEG_END ) {
                /* abnormal reading result; should not happen */
                while ( c != RI_ERR_EOF && !INCHI_INP_EOL(c) ) {
                    /* bypass to the end of line or file */
                    c = getInChIChar(pInp);
                }
                ret = (c == RI_ERR_EOF)? RI_ERR_EOF : RI_ERR_EOL; /* end of line */
                pLine->len = 0;
                pLine->c   = ret;
                break;
            }
            if ( c == RI_ERR_EOF ) {
                ret = RI_ERR_EOF; /* end of line */
                break;
            }
            if ( c == SEG_END ) {
                c = nGetInChISegment( pInp, pLine, szToken );
            }
            if ( c < 0 ) {
                goto exit_error; /* error */
            }
            if ( !pLine->len ) {
                ret = RI_ERR_EOL; /* end of line */
                break;
            }
            fst = UCINT pLine->str[0];
        }
        /* process the seqment */
        switch ( state ) {
        case AST_VERSION:
            /* Mobile H */
            bMobileH = TAUT_YES;
            ret = ParseAuxSegmentVersion( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_NUMBERS;
            break;
        case AST_MOBILE_H_NUMBERS:
            ret = ParseAuxSegmentNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = AST_MOBILE_H_ATOM_EQ;
            break;
        case AST_MOBILE_H_ATOM_EQ:
            ret = ParseAuxSegmentAtomEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_GROUP_EQ;
            break;
        case AST_MOBILE_H_GROUP_EQ:
            ret = ParseAuxSegmentGroupEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_SP3_INV;
            break;
        case AST_MOBILE_H_SP3_INV:
            ret = ParseAuxSegmentSp3Inv( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_SP3_INV_NUMBERS;
            break;
        case AST_MOBILE_H_SP3_INV_NUMBERS:
            ret = ParseAuxSegmentSp3InvNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_ISO_LAYER_FORK;
            break;
        case AST_MOBILE_H_ISO_LAYER_FORK:
            if ( !memcmp( pLine->str, "I:", 2 ) ) {
                state = AST_MOBILE_H_ISO_NUMBERS;
            } else
            if ( !memicmp( pLine->str, "F:", 2 ) ) { 
                state = AST_FIXED_H_NUMBERS;
                bMobileH = TAUT_NON;
            } else
            if ( /*bReconn == INCHI_BAS &&*/ !memicmp( pLine->str, "CRV:", 4 ) ) {
                state = AST_REVERSE_INFO_CRV;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "rA:", 3 ) ) {
                state = AST_REVERSE_INFO_ATOMS;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "R:", 3 ) ) {
                ret = 1;  /* read the next segment */
                state = AST_VERSION;
                bMobileH = TAUT_YES;
                bReconn = INCHI_REC;
            } else {
                ret = RI_ERR_SYNTAX;
            }
            break;
        /* Mobile H, isotopic */
        case AST_MOBILE_H_ISO_NUMBERS:
            ret = ParseAuxSegmentNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = AST_MOBILE_H_ISO_ATOM_EQ;
            break;
        case AST_MOBILE_H_ISO_ATOM_EQ:
            ret = ParseAuxSegmentAtomEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_ISO_GROUP_EQ;
            break;
        case AST_MOBILE_H_ISO_GROUP_EQ:
            ret = ParseAuxSegmentGroupEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_ISO_SP3_INV;
            break;
        case AST_MOBILE_H_ISO_SP3_INV:
            ret = ParseAuxSegmentSp3Inv( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_MOBILE_H_ISO_SP3_INV_NUMBERS;
            break;
        case AST_MOBILE_H_ISO_SP3_INV_NUMBERS:
            ret = ParseAuxSegmentSp3InvNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_LAYER_FORK;
            break;
        case AST_FIXED_H_LAYER_FORK:
            if ( !memicmp( pLine->str, "F:", 2 ) ) { 
                state = AST_FIXED_H_NUMBERS;
                bMobileH = TAUT_NON;
            } else
            if ( /*bReconn == INCHI_BAS &&*/ !memicmp( pLine->str, "CRV:", 4 ) ) {
                state = AST_REVERSE_INFO_CRV;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "rA:", 3 ) ) {
                state = AST_REVERSE_INFO_ATOMS;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "R:", 3 ) ) {
                ret = 1;  /* read the next segment */
                state = AST_VERSION;
                bMobileH = TAUT_YES;
                bReconn = INCHI_REC;
            } else {
                ret = RI_ERR_SYNTAX;
            }
            break;
        case AST_FIXED_H_NUMBERS:
            ret = ParseAuxSegmentNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = AST_FIXED_H_ATOM_EQ;
            break;
        case AST_FIXED_H_ATOM_EQ:
            ret = ParseAuxSegmentAtomEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_SP3_INV;
            break;
        case AST_FIXED_H_SP3_INV:
            ret = ParseAuxSegmentSp3Inv( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_SP3_INV_NUMBERS;
            break;
        case AST_FIXED_H_SP3_INV_NUMBERS:
            ret = ParseAuxSegmentSp3InvNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_ISO_LAYER_FORK;
            break;
        case AST_FIXED_H_ISO_LAYER_FORK:
            if ( !memcmp( pLine->str, "I:", 2 ) ) {
                state = AST_FIXED_H_ISO_NUMBERS;
            } else
            if ( /*bReconn == INCHI_BAS &&*/ !memicmp( pLine->str, "CRV:", 4 ) ) {
                state = AST_REVERSE_INFO_CRV;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "rA:", 3 ) ) {
                state = AST_REVERSE_INFO_ATOMS;
            } else
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "R:", 3 ) ) {
                ret = 1;  /* read the next segment */
                state = AST_VERSION;
                bMobileH = TAUT_YES;
                bReconn = INCHI_REC;
            } else {
                ret = RI_ERR_SYNTAX;
            }
            break;
        case AST_FIXED_H_ISO_NUMBERS:
            ret = ParseAuxSegmentNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = AST_FIXED_H_ISO_ATOM_EQ;
            break;
        case AST_FIXED_H_ISO_ATOM_EQ:
            ret = ParseAuxSegmentAtomEqu( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_SP3_INV;
            break;
        case AST_FIXED_H_ISO_SP3_INV:
            ret = ParseAuxSegmentSp3Inv( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_FIXED_H_ISO_SP3_INV_NUMBERS;
            break;
        case AST_FIXED_H_ISO_SP3_INV_NUMBERS:
            ret = ParseAuxSegmentSp3InvNumbers( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_REVERSE_INFO_CRV;
            break;
        case AST_REVERSE_INFO_CRV:
            ret = ParseAuxSegmentReverseCRV( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            /* state = (bReconn == INCHI_BAS)? AST_REVERSE_INFO_ATOMS : AST_RECONNECTED_LAYER_FORK;*/
            state = AST_REVERSE_INFO_ATOMS;
            break;
        case AST_REVERSE_INFO_ATOMS:
            ret = ParseAuxSegmentReverseAtoms( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_REVERSE_INFO_BONDS;
            break;
        case AST_REVERSE_INFO_BONDS:
            ret = ParseAuxSegmentReverseBonds( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_REVERSE_INFO_XYZ;
            break;
        case AST_REVERSE_INFO_XYZ:
            ret = ParseAuxSegmentReverseXYZ( pLine->str, bMobileH, &pXYZ, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = AST_RECONNECTED_LAYER_FORK;
            if ( ret > 0 ) {
                nLenXYZ = ret - 1;
            }
            break;
        case AST_RECONNECTED_LAYER_FORK:
            if ( bReconn == INCHI_BAS && !memicmp( pLine->str, "R:", 3 ) ) {
                ret = 1;  /* read the next segment */
                state = AST_VERSION;
                bMobileH = TAUT_YES;
                bReconn = INCHI_REC;
            } else {
                ret = RI_ERR_SYNTAX;
            }
            break;
        }
    } while( c >= 0 );
    
    ret = AddAuxSegmentCoord( ret, pXYZ, nLenXYZ, pInpInChI, nNumComponents );

exit_error:
    if ( pXYZ ) {
        inchi_free( pXYZ );
    }
    if ( ret >= 0 || c == RI_ERR_EOF || c == RI_ERR_EOL ) {
        pLine->len = 0;
    }

    return ret;
}
/******************************************************************************************************/
int ReadInChILine(INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, 
                  char **pStr, int *pState,
                  INChI *pInpInChI[INCHI_NUM][TAUT_NUM],
                  int nNumComponents[INCHI_NUM][TAUT_NUM],
                  REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM],
                  int s[INCHI_NUM][TAUT_NUM][2],
                  int *bStdFormat, 
                  int *bInputHasSaveOpt, unsigned char *inp_save_opt_bits)
{
    int   c, fst, ret=RI_ERR_ALLOC, len;
    int   bMobileH = TAUT_YES, bReconn = INCHI_BAS;
    const char szToken[] = INCHI_TOKEN;
    char *p;
    int  state=-1, prev_state=-1;
    int  bAbc = -1;  /* -1=> undefined, 0=> decimal, 1=> abc (compressed) */
    
    const int len_std_prefix=8;
    size_t k=0;
    unsigned char let1, let2;
    const char a2p[]="ABCDEFGHIJKLMNOP";

    /* memset( pLine, 0, sizeof( pLine[0] ) ); */
    *pState = 0;

next_line:
    INCHI_HEAPCHK
    /* Got "InChI=1/" */
    if ( pLine->len ) 
    {
        c = pLine->c;
    } 
    else 
    {
        INCHI_HEAPCHK
        c = nGetInChISegment( pInp, pLine, szToken );
        INCHI_HEAPCHK
    }
    if ( c == RI_ERR_EOF && !pLine->len && !pLine->str[0] ) 
    {
        ret = c;
        goto exit_function;
    }
    INCHI_HEAPCHK

    if ( pLine->len == 0 || (c != SEG_END && c != RI_ERR_EOF) || !(p = strstr(pLine->str, "InChI=1")) )
    {
        if ( pLine->str && pLine->str == strstr ( pLine->str, "Structure" ) ) 
        {
            if ( *pStr ) 
            {
                INCHI_HEAPCHK
                inchi_free( *pStr );
            }
            *pStr = pLine->str;
            /* bypass to the end of the 'Structure nnn' line */
            memset( pLine, 0, sizeof( pLine[0] ) );
            while ( c && !INCHI_INP_EOL(c) ) 
            {
                c = getInChIChar(pInp);
            }
            goto next_line;
        }
        /* bypass to the end of unrecognized line */
        while ( c != RI_ERR_EOF && !INCHI_INP_EOL(c) ) 
        {
            c = getInChIChar(pInp);
        }
        pLine->len = 0;
        INCHI_HEAPCHK
        goto next_line;
    }
    

    /* Check if got a standard InChI */	
    if ( ( pLine->len == len_std_prefix ) && (pLine->str[len_std_prefix-1]=='S') )
        *bStdFormat=1;
    else
        *bStdFormat=0;
    

    state=IST_MOBILE_H_FORMULA;
    ret  = 1; /* means read the next segment */
    do {
        /* read the next segment up to the '/' */
        INCHI_HEAPCHK
        if ( ret < 0 ) 
        {
            *pState = prev_state;
            break;
        }
        prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
        if ( 0 < ret ) 
        {
            /* read next segment */
            if ( c != RI_ERR_EOF && c != SEG_END ) 
            {
                /* abnormal reading result; should not happen */
                /* unless we got backslash-SaveOpt */
                if ( c=='\\' )
                {
                    /* May be SaveOpt */
                    *bInputHasSaveOpt = 1;
                }
                k = 0;
                while ( c != RI_ERR_EOF && !INCHI_INP_EOL(c) ) 
                {
                    /* bypass to the end of line or file */
                    c = getInChIChar(pInp);
                    k++;
                    if ( k==1 )       let1 = c;
                    else if ( k==2 )  let2 = c;
                }
                if ( k != 3) 
                {
                    /* not a valid SaveOpt which must be of two chars */
                    *bInputHasSaveOpt = 0;
                    let1 = let2 = '\0';
                }
                else        
                {
                    /* may be SaveOpt - analyze the content */
                    if  ( ( let2 >= 'A') && ( let2 <= 'D') )        /* letter-2 OK */
                    {
                        *bInputHasSaveOpt = 0;
                        *inp_save_opt_bits = 0;
                        for (k=0; k<16; k++)
                        {
                            if ( a2p[k] == let1)                    /* letter-1 OK */
                            {
                                *inp_save_opt_bits = (unsigned char) k;
                                *bInputHasSaveOpt = 1;
                                break;
                            }
                        }
                        if ( *bInputHasSaveOpt )
                        {
                            if ( let2=='B' || let2=='D' ) 
                                *inp_save_opt_bits |= SAVE_OPT_15T;
                            if ( let2=='C' || let2=='D' ) 
                                *inp_save_opt_bits |= SAVE_OPT_KET;
                        }
                    }
                }

                ret = (c == RI_ERR_EOF)? RI_ERR_EOF : RI_ERR_EOL; /* end of line */
                pLine->len = 0;
                pLine->c   = ret;
                break; /* exit */
            }
            if ( c == RI_ERR_EOF ) {
                ret = RI_ERR_EOF; /* end of line */
                break;
            }
            if ( c == SEG_END ) {
                c = nGetInChISegment( pInp, pLine, szToken );
            }
            if ( c < 0 ) {
                goto exit_error; /* error */
            }
            if ( !pLine->len ) {
                ret = RI_ERR_EOL; /* end of line */
                break;
            }
            fst = UCINT pLine->str[0];
        }
        /* process the seqment */
        switch ( state ) {

        /* Mobile H, M */                /* /  */
        case IST_MOBILE_H_FORMULA:
            bMobileH = TAUT_YES;
            ret = ParseSegmentFormula( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn] );
            state = IST_MOBILE_H_CONNECTIONS;
            break;
        case IST_MOBILE_H_CONNECTIONS:   /* /c */
            ret = ParseSegmentConnections( pLine->str, bMobileH, &pInpInChI[bReconn][bMobileH], &nNumComponents[bReconn][bMobileH], &bAbc );
            state = IST_MOBILE_H;
            break;
        case IST_MOBILE_H:               /* /h */
            ret = ParseSegmentMobileH( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], &bAbc );
            state = IST_MOBILE_H_CHARGE;
            break;
        case IST_MOBILE_H_CHARGE:        /* /q */
            ret = ParseSegmentCharge( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn] );
            state = IST_MOBILE_H_PROTONS;
            break;
        case IST_MOBILE_H_PROTONS:       /* /p */
            ret = ParseSegmentProtons( pLine->str, bMobileH, nNumProtons[bReconn], nNumComponents[bReconn] );
            state = IST_MOBILE_H_SP2;
            break;
        case IST_MOBILE_H_SP2:           /* /b */
            ret = ParseSegmentSp2( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_SP3;
            break;
        case IST_MOBILE_H_SP3:         /* t */
            ret = ParseSegmentSp3( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_SP3_M;
            break;
        case IST_MOBILE_H_SP3_M:       /* /m */
            ret = ParseSegmentSp3m( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = IST_MOBILE_H_SP3_S;
            break;
        case IST_MOBILE_H_SP3_S:       /* /s */ 
            ret = ParseSegmentSp3s( pLine->str, bMobileH, pInpInChI[bReconn], s[bReconn], nNumComponents[bReconn], state );
            state = IST_MOBILE_H_ISO_LAYER_FORK;
            break;
        case IST_MOBILE_H_ISO_LAYER_FORK:
            /* find layer type after M */
            ret = 0;
            switch( pLine->str[0] ) {
            case 'i':
                state = IST_MOBILE_H_ISO_ATOMS;  /* MI */
                break;
            case 'f':
                state = IST_FIXED_H_FORMULA; /* F */
                break;
            case 'r':
                state = IST_RECONNECTED_FORMULA; /* reconnected */
                break;
            default:
                ret = RI_ERR_SYNTAX;
            }
            if ( INCHI_INP_EOL(c) && ret == 0 && !pLine->str[1] ) {
                prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
                ret = RI_ERR_SYNTAX; /* empty layer /i or /f or /r at the end of InChI line */
            } else
            if ( !ret && state != IST_MOBILE_H_ISO_ATOMS ) {
                len = strlen( pLine->str );
                if ( len > 1 ) {
                    memmove( pLine->str, pLine->str+1, len );
                } else {
                    ret = 1; /* read the next segment */
                }
            }
            break;
        /* Mobile H, isotopic, MI */
        case IST_MOBILE_H_ISO_ATOMS:   /* i */
            ret = ParseSegmentIsoAtoms( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_ISO_EXCH_H;
            break;
        case IST_MOBILE_H_ISO_EXCH_H:  /* /i/h */
            ret = ParseSegmentIsoExchgH( pLine->str, bMobileH, nNumProtons[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_ISO_SP2;
            break;
        case IST_MOBILE_H_ISO_SP2:         /* /i/b */
            ret = ParseSegmentSp2( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_ISO_SP3;
            break;
        case IST_MOBILE_H_ISO_SP3:         /* /i/t */
            ret = ParseSegmentSp3( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_MOBILE_H_ISO_SP3_M;
            break;
        case IST_MOBILE_H_ISO_SP3_M:       /* /i/m */
            ret = ParseSegmentSp3m( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = IST_MOBILE_H_ISO_SP3_S;

            
            break;
        case IST_MOBILE_H_ISO_SP3_S:       /* /i/s */
            ret = ParseSegmentSp3s( pLine->str, bMobileH, pInpInChI[bReconn], s[bReconn], nNumComponents[bReconn], state );
            state = IST_FIXED_H_LAYER_FORK;
            break;
        case IST_FIXED_H_LAYER_FORK:
            /* find layer type after MI */
            ret = 0;
            switch( pLine->str[0] ) {
            case 'f':
                state = IST_FIXED_H_FORMULA; /* F */
                break;
            case 'r':
                state = IST_RECONNECTED_FORMULA; /* reconnected */
                break;
            default:
                ret = RI_ERR_SYNTAX;
            }
            if ( INCHI_INP_EOL(c) && ret == 0 && !pLine->str[1] ) {
                prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
                ret = RI_ERR_SYNTAX; /* empty layer /f or /r at the end of InChI line */
            } else
            if ( !ret ) {
                len = strlen( pLine->str );
                if ( len > 1 ) {
                    memmove( pLine->str, pLine->str+1, len );
                } else {
                    ret = 1; /* read the next segment */
                }
            }
            break;

        /* Fixed H, F */
        case IST_FIXED_H_FORMULA:
            bMobileH = TAUT_NON;
            ret = ParseSegmentFormula( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn] );
            state = IST_FIXED_H;
            break;
        case IST_FIXED_H:               /* /f/h */
            ret = ParseSegmentMobileH( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], &bAbc );
            state = IST_FIXED_H_CHARGE;
            break;
        case IST_FIXED_H_CHARGE:        /* /f/q */
            ret = ParseSegmentCharge( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn] );
            state = IST_FIXED_H_SP2;
            break;
        case IST_FIXED_H_SP2:           /* /f/b */
            ret = ParseSegmentSp2( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_SP3;
            break;
        case IST_FIXED_H_SP3:         /* /f/t */
            ret = ParseSegmentSp3( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_SP3_M;
            break;
        case IST_FIXED_H_SP3_M:       /* /f/m */
            ret = ParseSegmentSp3m( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = IST_FIXED_H_SP3_S;
            break;
        case IST_FIXED_H_SP3_S:       /* /f/s */
            ret = ParseSegmentSp3s( pLine->str, bMobileH, pInpInChI[bReconn], s[bReconn], nNumComponents[bReconn], state );
            state = IST_FIXED_H_PERMUTATION;
            break;
        case IST_FIXED_H_PERMUTATION:  /* /f/o */
            ret = ParseSegmentPerm( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_ISO_LAYER_FORK;
            break;
        case IST_FIXED_H_ISO_LAYER_FORK:
            /* find layer type after M */
            ret = 0;
            switch( pLine->str[0] ) {
            case 'i':
                state = IST_FIXED_H_ISO_ATOMS;  /* FI */
                break;
            case 'r':
                state = IST_RECONNECTED_FORMULA; /* reconnected */
                break;
            default:
                ret = RI_ERR_SYNTAX;
            }
            if ( INCHI_INP_EOL(c) && ret == 0 && !pLine->str[1] ) {
                prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
                ret = RI_ERR_SYNTAX; /* empty layer /i or /r at the end of InChI line */
            } else
            if ( !ret && state != IST_FIXED_H_ISO_ATOMS ) {
                len = strlen( pLine->str );
                if ( len > 1 ) {
                    memmove( pLine->str, pLine->str+1, len );
                } else {
                    ret = 1; /* read the next segment */
                }
            }
            break;

        /* Fixed H, isotopic, FI */
        case IST_FIXED_H_ISO_ATOMS:   /* /f/i */
            ret = ParseSegmentIsoAtoms( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_ISO_SP2;
            break;
        case IST_FIXED_H_ISO_SP2:         /* /f/i/b */
            ret = ParseSegmentSp2( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_ISO_SP3;
            break;
        case IST_FIXED_H_ISO_SP3:         /* /f/i/t */
            ret = ParseSegmentSp3( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_FIXED_H_ISO_SP3_M;
            break;
        case IST_FIXED_H_ISO_SP3_M:       /* /f/i/m */
            ret = ParseSegmentSp3m( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state );
            state = IST_FIXED_H_ISO_SP3_S;
            break;
        case IST_FIXED_H_ISO_SP3_S:       /* /f/i/s */
            ret = ParseSegmentSp3s( pLine->str, bMobileH, pInpInChI[bReconn], s[bReconn], nNumComponents[bReconn], state );
            state = IST_FIXED_H_ISO_PERMUTATION;
            break;
        case IST_FIXED_H_ISO_PERMUTATION:  /* /f/i/o */
            ret = ParseSegmentPerm( pLine->str, bMobileH, pInpInChI[bReconn], nNumComponents[bReconn], state, &bAbc );
            state = IST_RECONNECTED_LAYER_FORK;
            break;
        case IST_RECONNECTED_LAYER_FORK:
            /* find layer type after FI */
            ret = 0;
            switch( pLine->str[0] ) {
            case 'r':
                state = IST_RECONNECTED_FORMULA; /* reconnected */
                break;
            default:
                ret = RI_ERR_SYNTAX;
            }
            if ( INCHI_INP_EOL(c) && ret == 0 && !pLine->str[1] ) {
                prev_state = state + (bReconn? IST_HAPPENED_IN_RECMET : 0);
                ret = RI_ERR_SYNTAX; /* empty layer /r at the end of InChI line */
            } else
            if ( !ret ) {
                len = strlen( pLine->str );
                if ( len > 1 ) {
                    memmove( pLine->str, pLine->str+1, len );
                } else {
                    ret = 1; /* read the next segment */
                }
            }
            break;
        case IST_RECONNECTED_FORMULA:
            bReconn = INCHI_REC;
            bMobileH = TAUT_YES;
            state   = IST_MOBILE_H_FORMULA;
            break;
        }


    } while( c >= 0 );

exit_function:;
exit_error:;

    INCHI_HEAPCHK

    if ( ret >= 0 || c == RI_ERR_EOF || c == RI_ERR_EOL ) {
        pLine->len = 0;
    }
    return ret;
}
/****************************************************************************************/
int ParseSegmentIsoExchgH( const char *str, int bMobileH, REM_PROTONS nNumProtons[], int pnNumComponents[], int state, int *pbAbc )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    const char *p, *q, *pStart, *pEnd;
    int  ret=0, num, i, i_prev;
    static char abc_h[] = "hdt";

    if ( str[0] != 'h' )
        return 0;

    pStart = str+1;

    if ( !(bMobileH==TAUT_YES && state == IST_MOBILE_H_ISO_EXCH_H ) ) {
        return RI_ERR_PROGR; /* program error */
    }

    if ( !(pEnd = strchr( pStart, ';' )) ) {
        pEnd = pStart + strlen(pStart);
    } else {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    p = pStart;

    if ( p < pEnd && *pbAbc == -1 ) {
        /* check if compressed InChI */
        /* compressed:    /hNtNdNh where N is a decimal number */
        /* uncompressed:  /hT[n]D[n]H[n] where n > 1 is a decimal number */ 
        *pbAbc = isdigit( UCINT *p)? 1 : 0;
    }

    if ( *pbAbc == 1 ) {
        i_prev = (int)sizeof(abc_h);
        while ( p < pEnd ) {
            num = (int)inchi_strtol( p, &q, 10 );
            if ( 0 >= num  || p == q || q >= pEnd ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            p = strchr( abc_h, *q);
            if ( p && (i=p-abc_h) < i_prev ) {
                nNumProtons[bMobileH].nNumRemovedIsotopicH[i] = (NUM_H)num;
                p = q+1;
                i_prev = i;
            } else {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
        }
    } else {
        if ( *p == 'T' ) {
            nNumProtons[bMobileH].nNumRemovedIsotopicH[2] = 1;
            p ++;
            if ( isdigit( UCINT p[0]) ) {
                nNumProtons[bMobileH].nNumRemovedIsotopicH[2] = (NUM_H)inchi_strtol( p, &q, 10 );
                p = q;
            }
        }
        if ( *p == 'D' ) {
            nNumProtons[bMobileH].nNumRemovedIsotopicH[1] = 1;
            p ++;
            if ( isdigit( UCINT p[0]) ) {
                nNumProtons[bMobileH].nNumRemovedIsotopicH[1] = (NUM_H)inchi_strtol( p, &q, 10 );
                p = q;
            }
        }
        if ( *p == 'H' ) {
            nNumProtons[bMobileH].nNumRemovedIsotopicH[0] = 1;
            p ++;
            if ( isdigit( UCINT p[0]) ) {
                nNumProtons[bMobileH].nNumRemovedIsotopicH[0] = (NUM_H)inchi_strtol( p, &q, 10 );
                p = q;
            }
        }
    }
    if ( p != pEnd ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentPerm( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state, int *pbAbc )
{
    int nNumComponents, iComponent1, iComponent2, numTrans;
    const char *p, *q, *pStart, *pEnd, *pPermStart, *pPermEnd;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH]; /* bMobileH should be TAUT_NON = 0 */
    INChI tmp;
    int   base = 10;


    if ( str[0] != 'o' )
        return 0;

    pStart = str+1;
    nNumComponents = ppnNumComponents[bMobileH];

    if ( !(bMobileH==TAUT_NON && ( state == IST_FIXED_H_PERMUTATION || state == IST_FIXED_H_ISO_PERMUTATION) ) ) {
        return RI_ERR_PROGR; /* program error */
    }

    if ( !(pEnd = strchr( pStart, ';' )) ) {
        pEnd = pStart + strlen(pStart);
    } else {
        return RI_ERR_SYNTAX; /* syntax error */
    }
    while( pStart < pEnd ) {
        /* cycle over components; rearrange Fixed H components in order of Mobile H components */
        /* if /o(1,2,3) then reaarange Fixed H components in this way: tmp<-1, 1<-2, 2<-3, 3<-tmp */
        if ( *pStart != '(' ) {
            ret = RI_ERR_SYNTAX;
            goto exit_function;
        }
        pPermStart = pStart + 1;
        memset( &tmp, 0, sizeof(tmp) );  /* initialization 2006-03 */
        if ( !(pPermEnd   = strchr( pPermStart, ')' )) || pPermEnd == pPermStart ) {
            ret = RI_ERR_SYNTAX;
            goto exit_function;
        }
        
        if ( pPermStart < pPermEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = isupper( UCINT *pPermStart)? 1 : 0;
        }
        base = (*pbAbc==1)? ALPHA_BASE : 10;

        /* permutation cycle */
        if ( *pbAbc == 1 ) {
            for ( p = pPermStart, iComponent2 = numTrans = 0; p < pPermEnd; iComponent2 = iComponent1, p = q ) {
                /* get first atom number */
                if ( 0 >= (iComponent1 = (int)inchi_strtol( p, &q, base )) || iComponent1 > nNumComponents ) {
                    ret = RI_ERR_SYNTAX;  /* syntax error */
                    goto exit_function;
                }
                if ( iComponent2 ) {
                    pInChI[iComponent2-1] = pInChI[iComponent1-1];
                    numTrans ++;
                } else {
                    tmp = pInChI[iComponent1-1]; /* on the 1st pass save Component1 */
                }
            }
        } else {
            for ( p = pPermStart, iComponent2 = numTrans = 0; p < pPermEnd; iComponent2 = iComponent1, p = q + (*q==',') ) {
                /* get first atom number */
                if ( !isdigit( UCINT *p ) ) {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                if ( !(iComponent1 = (int)inchi_strtol( p, &q, 10 )) || iComponent1 > nNumComponents ) {
                    ret = RI_ERR_SYNTAX;  /* syntax error */
                    goto exit_function;
                }
                if ( iComponent2 ) {
                    pInChI[iComponent2-1] = pInChI[iComponent1-1];
                    numTrans ++;
                } else {
                    tmp = pInChI[iComponent1-1]; /* on the 1st pass save Component1 */
                }
            }
        }
        pInChI[iComponent2-1] = tmp;
        if ( !numTrans || p != pPermEnd ) {
            ret = RI_ERR_SYNTAX;
            goto exit_function;
        } else {
            pStart = p+1;
        }
    }
    ret = 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentIsoAtoms( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state, int *pbAbc )
{
    int i, mpy_component, val;
    int nNumComponents, iComponent, len, iAtom;
    AT_NUMB nAtom1;
    const char *p, *q, *t, *pStart, *pEnd, *r;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH];
    INChI *pInChIFrom=NULL;
    INChI_IsotopicAtom **pIsotopicAtom = NULL;
    INChI_IsotopicAtom isoAtom;

    const char   mult_type[]       = "mnMNe";
    const char   parity_type[]     = "-+TDH";
    int    bIsoFrom, nCpyType = CPY_ISO_AT;
    int    base = 10;

    if ( str[0] != 'i' )
        return 0;

    pStart = str+1;
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    if ( !((bMobileH==TAUT_YES && state == IST_MOBILE_H_ISO_ATOMS) ||
           (bMobileH==TAUT_NON && state == IST_FIXED_H_ISO_ATOMS) ) ) {
        return RI_ERR_PROGR; /* program error */
    }
    if ( !*pStart ) {
        return nNumComponents+1; /* no isotopic atoms */
    }

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
            if ( p != q ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#if (FIX_DALKE_BUGS == 1)
            if ( iComponent + mpy_component > nNumComponents ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#endif
            p ++; /* move to the 1st character of the component */
        } else
        if ( ((isdigit(*pStart) &&
             0 < (val = (int)inchi_strtol( pStart, &q, 10))) ||
             (q = pStart, val=1))&&
             (t=strchr(mult_type, *q)) && q+1 == pEnd ) {
            /* process the abbreviation */
            ret = 0;
#if (FIX_DALKE_BUGS == 1)
            if ( iComponent + val > nNumComponents ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#endif
            bIsoFrom   = 0;
            switch( bMobileH ) {
            case TAUT_YES:
                ret = RI_ERR_SYNTAX;
                break;
            case TAUT_NON:
                if ( *q == 'm' ) {
                    /* copy from mobile H to fixed H */
                    pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                } else
                if ( *q == 'e' ) {
                    /* copy from mobile H to isotopic mobile H */
                    pInChIFrom = pInChI;
                    bIsoFrom   = -1; /* empty */
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                }
                break;
            default:
                ret = RI_ERR_SYNTAX; 
                break;
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
            /* copy */
            for ( i = 0; i < val; i ++ ) {
                ret = CopySegment( pInChI+iComponent+i, pInChIFrom+iComponent+i, nCpyType, 0, bIsoFrom );
                if ( !ret ) {
                    ret = RI_ERR_SYNTAX;
                }
                if ( ret < 0 ) {
                    goto exit_function;
                }
            }
            iComponent += val;
            /* continue to the next component(s) */
            if ( *pEnd ) {
                pStart = pEnd+1;
                continue;
            } else {
                break;
            }
        } else {
            mpy_component = 1;
            p = pStart;
        }
        pStart = p;
        pIsotopicAtom = &pInChI[iComponent].IsotopicAtom;
        if ( *pIsotopicAtom ) {
            ret = RI_ERR_PROGR; /* program error */
            goto exit_function;
        }

        if ( p < pEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = isupper( UCINT *p)? 1 : 0;
        }
        base = (*pbAbc==1)? ALPHA_BASE : 10;


one_more_time:
        if ( *pbAbc == 1 ) {
            /* process the componnt: At[+/-Charge]TDH,... */
            /* pass 1: find number of stereoatoms */
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++ ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &p, base );
                if ( !nAtom1 ||
                     nAtom1 > pInChI[iComponent].nNumberOfAtoms ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                memset( &isoAtom, 0, sizeof(isoAtom) );
                isoAtom.nAtomNumber = nAtom1;
                isoAtom.nIsoDifference = (NUM_H)inchi_strtol( p, &q, 10 ); /* alway in abc */
                if ( p == q ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q;
                if ( *p == 't' ) {
                    isoAtom.nNum_T = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_T = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( *p == 'd' ) {
                    isoAtom.nNum_D = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_D = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( *p == 'h' ) {
                    isoAtom.nNum_H = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_H = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( p > pEnd || (!isoAtom.nIsoDifference && !isoAtom.nNum_T && !isoAtom.nNum_D && !isoAtom.nNum_H) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( *pIsotopicAtom ) {
                    pIsotopicAtom[0][iAtom] = isoAtom;
                }
            }
        } else {
            /* process the componnt: At[+/-Charge]TDH,... */
            /* pass 1: find number of stereoatoms */
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++ ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                p = q;
                if ( !nAtom1 ||
                     nAtom1 > pInChI[iComponent].nNumberOfAtoms ||
                     !(r = strchr( parity_type, *p) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                memset( &isoAtom, 0, sizeof(isoAtom) );
                isoAtom.nAtomNumber = nAtom1;
                if ( p[0] == '+' && isdigit( UCINT p[1]) ) {
                    isoAtom.nIsoDifference = (NUM_H)inchi_strtol( p+1, &q, 10 );
                    if ( isoAtom.nIsoDifference >= 0 ) isoAtom.nIsoDifference ++;
                    p = q;
                } else
                if ( p[0] == '-' && isdigit( UCINT p[1]) ) {
                    isoAtom.nIsoDifference = -(NUM_H)inchi_strtol( p+1, &q, 10 );
                    if ( isoAtom.nIsoDifference == 0 ) isoAtom.nIsoDifference ++;
                    p = q;
                }
                if ( *p == 'T' ) {
                    isoAtom.nNum_T = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_T = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( *p == 'D' ) {
                    isoAtom.nNum_D = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_D = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( *p == 'H' ) {
                    isoAtom.nNum_H = 1;
                    p ++;
                    if ( isdigit( UCINT *p) ) {
                        isoAtom.nNum_H = (NUM_H)inchi_strtol( p, &q, 10 );
                        p = q;
                    }
                }
                if ( !isoAtom.nIsoDifference && !isoAtom.nNum_T && !isoAtom.nNum_D && !isoAtom.nNum_H ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( p < pEnd ) {
                    if ( *p == ',' ) {
                        p ++;
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                }
                if ( *pIsotopicAtom ) {
                    pIsotopicAtom[0][iAtom] = isoAtom;
                }
            }
        }
        if ( p != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }

        if ( !*pIsotopicAtom ) {
            /* end of the 1st pass */
            len = iAtom;
            /* memory allocation */
            if ( !(*pIsotopicAtom = (INChI_IsotopicAtom *) inchi_calloc( len+1, sizeof(**pIsotopicAtom) ) ) ) {
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
            goto one_more_time; /* goto the 2nd pass */
        } else {
            /* 2nd pass */
            if ( len != iAtom ) {
                ret = RI_ERR_PROGR; /* program error */
                goto exit_function;
            }
            pInChI[iComponent].nNumberOfIsotopicAtoms = len;
        }

        /* multiplier */
        for ( i = 1; i < mpy_component; i ++ ) {
            ret = CopySegment( pInChI+iComponent+i, pInChI+iComponent, nCpyType, 0, 0 );
            if ( !ret ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
        }

        iComponent += mpy_component;
        if ( *pEnd ) {
            pStart = pEnd+1;
            continue;
        } else {
            break;
        }
    }
    if ( nNumComponents != iComponent ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentSp3s( const char *str, int bMobileH, INChI *pInpInChI[], int s[TAUT_NUM][2], int ppnNumComponents[], int state )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int nNumComponents, iComponent, val;
    const char *p, *q, *pStart, *pEnd;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH];
    INChI_Stereo **pStereo = NULL;

    int   bIso = (state==IST_MOBILE_H_ISO_SP3_S || state==IST_FIXED_H_ISO_SP3_S);
    
    if ( !bIso && state != IST_MOBILE_H_SP3_S && state != IST_FIXED_H_SP3_S ) {
        return RI_ERR_PROGR; /* program error */
    }

    if ( str[0] != 's' )
        return 0;

    pStart = str+1;
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    /*if ( !(pEnd = strchr( pStart, ';' )) )*/ /* 2007-09-25 DT */
    if ( !(pEnd = strchr( pStart, '/' )) ){
        pEnd = pStart + strlen(pStart);
    } else {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    p = pStart;
    if ( pEnd == pStart ) {
        /* create empty sp3 segment */
        int len = 0;
        s[bMobileH][bIso] = NO_VALUE_INT; /* empty */
        /* create empty sp3 segment */
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            pStereo = bIso? &pInChI[iComponent].StereoIsotopic : &pInChI[iComponent].Stereo;
            if ( !*pStereo ) {
                if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                    ret = RI_ERR_ALLOC; /* memory allocation failed */
                    goto exit_function;
                }
            }
            pStereo[0]->nCompInv2Abs = 0;  /* deliberately empty */

            if ( pStereo[0]->nNumberOfStereoCenters ) {
                ret = RI_ERR_SYNTAX; /* syntax error: "/s" without a digit describes "no stereo" */
                goto exit_function;
            }
            /* allocate empty sp3 stereo */
            if ( (!pStereo[0]->t_parity &&
                 !(pStereo[0]->t_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pStereo[0]->b_parity[0]) ) )) ||
                 (!pStereo[0]->nNumber &&
                 !(pStereo[0]->nNumber = (AT_NUMB *)inchi_calloc( len+1, sizeof(pStereo[0]->nNumber[0]) ) )) ) {
                /* cleanup */
                if ( pStereo[0]->t_parity ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->t_parity );
                    pStereo[0]->t_parity = NULL;
                }
                if ( pStereo[0]->nNumber ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->nNumber );
                    pStereo[0]->nNumber = NULL;
                }
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
        }
        ret = nNumComponents+1;
    } else {
        val = (int)inchi_strtol( p, &q, 10 );
        if ( q == pEnd && 1 <= val && val <= 3 ) {
            s[bMobileH][bIso] = val;
            ret = nNumComponents+1;
        } else {
            ret = RI_ERR_SYNTAX; /* syntax error */
        }
    }
exit_function:
    return ret;
}
/****************************************************************************************/
int bIsSp3LayerNotEmpty( INChI *pInpInChI[], int bMobileH, int bIso, int nNumComponents )
{
    INChI        *pInChI;
    INChI_Stereo *pStereo;
    int           iComponent, num_not_empty = 0;

    if ( pInpInChI[bMobileH] ) {
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            pInChI  = pInpInChI[bMobileH] + iComponent;
            if ( pInChI->bDeleted || !pInChI->nNumberOfAtoms ) {
                continue;
            }
            pStereo = bIso? pInChI->StereoIsotopic : pInChI->Stereo;
            if ( pStereo && pStereo->nNumberOfStereoCenters > 0 && pStereo->nNumber && pStereo->t_parity ) {
                num_not_empty ++;
            }
        }
    }
    return num_not_empty;
}
/****************************************************************************************/
int ParseSegmentSp3m( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int nNumComponents, iComponent;
    const char *p, *pStart, *pEnd;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH];
    INChI_Stereo **pStereo = NULL;

    int   bIso = (state==IST_MOBILE_H_ISO_SP3_M || state==IST_FIXED_H_ISO_SP3_M);
    
    if ( !bIso && state != IST_MOBILE_H_SP3_M && state != IST_FIXED_H_SP3_M ) {
        return RI_ERR_PROGR; /* program error */
    }
    nNumComponents = ppnNumComponents[bMobileH];

    if ( str[0] != 'm' ) {
        /* /m is missing: check whether we have to inherit /m from a preceding stereo layer */
        INChI_Stereo *pStereoFrom, *pStereoTo;
        INChI        *pInChIFrom;
        int          nNumCopied = 0, bMobileHFrom=-1, bIsoFrom=-1;
        if ( bMobileH && !bIso ) {
            return 0; /* Main non-isotopic cannot inherit: it has no preceding layer */
        } else
        if ( !bMobileH && !bIso ) {
            /* fixed-H non-isotopic (F) inherits from Mobile-H non-isotopic (M) */
            bMobileHFrom = TAUT_YES;
            bIsoFrom     = 0;
        } else
        if ( bMobileH && bIso ) {
            /* Mobile-H isotopic (MI) inherits from Mobile-H non-isotopic (M) */
            bMobileHFrom = TAUT_YES;
            bIsoFrom     = 0;
        } else
        if ( !bMobileH && bIso ) {
            /* Fixed-H isotopic (FI) inherits from Fixed-H non-isotopic (F) */
            bMobileHFrom = TAUT_NON;
            bIsoFrom     = 0;
            /* if Sp3 is empty in F as well as in M, then inherit from MI */
            if ( !bIsSp3LayerNotEmpty( pInpInChI, TAUT_NON, 0, ppnNumComponents[TAUT_NON /*bMobileH*/] ) /* F */ &&
                 !bIsSp3LayerNotEmpty( pInpInChI, TAUT_YES, 0, ppnNumComponents[TAUT_YES /*bMobileH*/] ) /* M */ ) {
                bMobileHFrom = TAUT_YES;
                bIsoFrom     = 1;
            }
        }
        if ( bMobileHFrom < 0 || bIsoFrom < 0 ) {
            return RI_ERR_PROGR;
        }
        if ( !bIsSp3LayerNotEmpty( pInpInChI, bMobileHFrom, bIsoFrom, ppnNumComponents[/*bMobileH*/ bMobileHFrom] ) ) {
            /* nothing to copy; check whether it should have inherited from a preceding layer */
            if ( (!bMobileHFrom && bIsoFrom) || (bMobileHFrom && !bIsoFrom) ) {
                /* MI or F inherit stereo from M */
                bMobileHFrom = TAUT_YES;
                bIsoFrom     = 0;
                if ( !bIsSp3LayerNotEmpty( pInpInChI, bMobileHFrom, bIsoFrom, ppnNumComponents[bMobileHFrom /*bMobileH*/] ) ) {
                    return 0;
                }
            } else {
                return 0;
            }
        }
        nNumComponents = inchi_min( ppnNumComponents[bMobileH], ppnNumComponents[bMobileHFrom] );
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            pInChIFrom = pInpInChI[bMobileHFrom] + iComponent;
            pInChI     = pInpInChI[bMobileH]     + iComponent;
            if ( pInChIFrom->nNumberOfAtoms > 0 && !pInChIFrom->bDeleted &&
                 pInChI->nNumberOfAtoms > 0    && !pInChI->bDeleted        ) {
                pStereoFrom = bIsoFrom? pInChIFrom->StereoIsotopic : pInChIFrom->Stereo;
                pStereoTo    = bIso?     pInChI->StereoIsotopic     : pInChI->Stereo;
                if ( pStereoFrom && pStereoTo ) {
                    pStereoTo->nCompInv2Abs = pStereoFrom->nCompInv2Abs;
                    nNumCopied ++;
                }
            }
        }
        return 0; /* return value > 0 means the non-/m segment has been processed here */
    }

    pStart = str+1;
    iComponent = 0;

    /*if ( !(pEnd = strchr( pStart, ';' )) )*/ /* 2007-09-25 DT */
    if ( !(pEnd = strchr( pStart, '/' )) ) {
        pEnd = pStart + strlen(pStart);
    } else {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    p = pStart;
    if ( pEnd == pStart ) {
        /* create empty sp3 segment */
        int len = 0;
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            INChI *pIsoInChI = &pInChI[iComponent];
            pStereo = bIso? &pIsoInChI->StereoIsotopic : &pIsoInChI->Stereo;
            if ( !*pStereo ) {
                if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                    ret = RI_ERR_ALLOC; /* memory allocation failed */
                    goto exit_function;
                }
            }
            pStereo[0]->nCompInv2Abs = NO_VALUE_INT;  /* deliberately empty */
#ifdef NEVER            
            if ( pStereo[0]->nNumberOfStereoCenters ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#endif
            /* allocate empty sp3 stereo */
            if ( (!pStereo[0]->t_parity &&
                 !(pStereo[0]->t_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pStereo[0]->b_parity[0]) ) )) ||
                 (!pStereo[0]->nNumber &&
                 !(pStereo[0]->nNumber = (AT_NUMB *)inchi_calloc( len+1, sizeof(pStereo[0]->nNumber[0]) ) )) ) {
                /* cleanup */
                if ( pStereo[0]->t_parity ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->t_parity );
                    pStereo[0]->t_parity = NULL;
                }
                if ( pStereo[0]->nNumber ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->nNumber );
                    pStereo[0]->nNumber = NULL;
                }
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
        }
        ret = nNumComponents+1;
    } else {
        while( p < pEnd && iComponent < nNumComponents ) {
            /* cycle over components */
            pStereo = bIso? &pInChI[iComponent].StereoIsotopic : &pInChI[iComponent].Stereo;
            if ( *p != '.' && !*pStereo ) {
                if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                    ret = RI_ERR_ALLOC; /* memory allocation failed */
                    goto exit_function;
                }
            }
            switch( *p ) {
            case '1':
                pStereo[0]->nCompInv2Abs = -1;
                break;
            case '0':
                pStereo[0]->nCompInv2Abs = 1;
                break;
            case '.':
                if ( *pStereo ) {
                    pStereo[0]->nCompInv2Abs = 0;
                }
                break;
            default:
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            iComponent ++;
            p ++;
        }
        if ( p != pEnd || iComponent != nNumComponents ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        ret = nNumComponents+1;
    }
exit_function:
    return ret;
}
/****************************************************************************************/
int ParseSegmentSp3( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state, int *pbAbc )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int i, mpy_component, val;
    int nNumComponents, iComponent, len, iAtom;
    AT_NUMB nAtom1;
    int     atomParity;
    const char *p, *q, *t, *pStart, *pEnd, *r;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH];
    INChI *pInChIFrom=NULL;
    /*
    INChI_Stereo *Stereo = NULL;
    INChI_Stereo *StereoOther = NULL;
    */
    INChI_Stereo **pStereo = NULL;

    const char   mult_type[]   = "mnMNe";
    const char   parity_type[] = "-+u?";
    int   bIsoTo, bIsoFrom, nCpyType = CPY_SP3;
    int   bIso = (state==IST_MOBILE_H_ISO_SP3 || state==IST_FIXED_H_ISO_SP3);
    int   base = 10;
    
    if ( !bIso && state != IST_MOBILE_H_SP3 && state != IST_FIXED_H_SP3 ) {
        return RI_ERR_PROGR; /* program error */
    }

    if ( str[0] != 't' )
        return 0;

    pStart = str+1;
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    if ( !*pStart ) {
        /* create empty sp3 segment */
        int len0 = 0;
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            INChI *pIsoInChI = &pInChI[iComponent];
            pStereo = bIso? &pIsoInChI->StereoIsotopic : &pIsoInChI->Stereo;
            if ( !*pStereo ) {
                if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                    ret = RI_ERR_ALLOC; /* memory allocation failed */
                    goto exit_function;
                }
            }
            /* allocate empty sp3 stereo */
            if ( (!pStereo[0]->b_parity &&
                 !(pStereo[0]->b_parity = (S_CHAR *)inchi_calloc( len0+1, sizeof(pStereo[0]->b_parity[0]) ) )) ||
                 (!pStereo[0]->nBondAtom1 &&
                 !(pStereo[0]->nBondAtom1 = (AT_NUMB *)inchi_calloc( len0+1, sizeof(pStereo[0]->nBondAtom1[0]) ) )) ||
                 (!pStereo[0]->nBondAtom2 &&
                 !(pStereo[0]->nBondAtom2 = (AT_NUMB *)inchi_calloc( len0+1, sizeof(pStereo[0]->nBondAtom2[0]) ) )) ) {
                /* cleanup */
                if ( pStereo[0]->b_parity ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->b_parity );
                    pStereo[0]->b_parity = NULL;
                }
                if ( pStereo[0]->nBondAtom1 ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->nBondAtom1 );
                    pStereo[0]->nBondAtom1 = NULL;
                }
                if ( pStereo[0]->nBondAtom2 ) {
                    INCHI_HEAPCHK
                    inchi_free( pStereo[0]->nBondAtom2 );
                    pStereo[0]->nBondAtom2 = NULL;
                }
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
            pStereo[0]->nCompInv2Abs = NO_VALUE_INT;
        }
        ret = nNumComponents+1;
        goto exit_function;
    }

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        if ( ((isdigit(*pStart) &&
             0 < (val = (int)inchi_strtol( pStart, &q, 10))) ||
             (q = pStart, val=1))&&
             (t=strchr(mult_type, *q)) && q+1 == pEnd ) {
            /* process the abbreviation */
            ret = 0;
#if (FIX_DALKE_BUGS == 1)
        if ( iComponent + val > nNumComponents ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
#endif
            switch( bMobileH ) {
            case TAUT_YES:
                switch( state ) {
                case IST_MOBILE_H_ISO_SP3:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'e' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = -1; /* empty */
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                default:
                    ret = RI_ERR_SYNTAX; 
                    break;
                }
                break;
            case TAUT_NON:
                switch( state ) {
                case IST_FIXED_H_SP3:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to fixed H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 0;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'e' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = -1; /* empty */
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                case IST_FIXED_H_ISO_SP3:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to fixed isotopic H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'M' ) {
                        /* copy from isotopic mobile H to fixed isotopic H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 1;
                        bIsoFrom   = 1;
                    } else
                    if ( *q == 'n' ) {
                        /* copy from fixed H to fixed isotopic H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'e' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = -1; /* empty */
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                default:
                    ret = RI_ERR_SYNTAX; 
                    break;
                }
                break;

            default:
                ret = RI_ERR_SYNTAX; 
                break;

            }
            if ( ret < 0 ) {
                goto exit_function;
            }
            /* copy */
            for ( i = 0; i < val; i ++ ) {
                ret = CopySegment( pInChI+iComponent+i, pInChIFrom+iComponent+i, nCpyType, bIsoTo, bIsoFrom );
                if ( !ret ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                }
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( bIsoFrom >= 0 ) {
                    INChI_Stereo *pStereoTo = bIsoTo? pInChI[iComponent+i].StereoIsotopic : pInChI[iComponent+i].Stereo;
                    if ( pStereoTo ) {
                        pStereoTo->nCompInv2Abs = NO_VALUE_INT; /* in case there in no /m segment after this */
                    }
                }
            }
            
            mpy_component = val;
            goto end_main_cycle;

        } else
        /* regular multiplier */
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
            if ( p != q ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            p ++; /* move to the 1st character of the component */
        } else {
            mpy_component = 1;
            p = pStart;
        }
#if (FIX_DALKE_BUGS == 1)
        if ( iComponent + mpy_component > nNumComponents ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
#endif
        pStart = p;
        if ( p < pEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = isupper( UCINT *p)? 1 : 0;
        }
        base = (*pbAbc==1)? ALPHA_BASE : 10;
        /* process the componnt: at1p,at1p,... */
        /* pass 1: find number of stereoatoms */
        if ( *pbAbc == 1 ) {
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++ ) {
                if ( (nAtom1 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (atomParity = (int)inchi_strtol( p, &p, 10),
                      AB_MIN_KNOWN_PARITY <= atomParity && atomParity <= AB_MAX_KNOWN_PARITY) ) {
                    ; /* okay */
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( nAtom1 > pInChI[iComponent].nNumberOfAtoms ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        } else {
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++, p += (*p == ',') ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                p = q+1;
                if ( !nAtom1 ||
                     nAtom1 > pInChI[iComponent].nNumberOfAtoms ||
                     !(r = strchr( parity_type, *q) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        }
        if ( p != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        len = iAtom;

        /* memory allocation */

        pStereo = bIso? &pInChI[iComponent].StereoIsotopic : &pInChI[iComponent].Stereo;

        if ( !*pStereo ) {
            if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
        }
        if ( pStereo[0]->t_parity || pStereo[0]->nNumberOfStereoCenters ||
             pStereo[0]->nNumber ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        /* allocate sp3 stereo */
        if ( !(pStereo[0]->t_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pStereo[0]->b_parity[0]) ) ) ||
             !(pStereo[0]->nNumber = (AT_NUMB *)inchi_calloc( len+1, sizeof(pStereo[0]->nNumber[0]) ) ) ) {
            /* cleanup */
            if ( pStereo[0]->t_parity ) {
                INCHI_HEAPCHK
                inchi_free( pStereo[0]->t_parity );
                pStereo[0]->t_parity = NULL;
            }
            if ( pStereo[0]->nNumber ) {
                INCHI_HEAPCHK
                inchi_free( pStereo[0]->nNumber );
                pStereo[0]->nNumber = NULL;
            }
            ret = RI_ERR_ALLOC; /* memory allocation failed */
            goto exit_function;
        }

        /* pass 2: store stereocenters */
        if ( *pbAbc == 1 ) {
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++ ) {
                if ( (nAtom1 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (atomParity = (int)inchi_strtol( p, &p, 10),
                      AB_MIN_KNOWN_PARITY <= atomParity && atomParity <= AB_MAX_KNOWN_PARITY) ) {
                    ; /* okay */
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( nAtom1 > pInChI[iComponent].nNumberOfAtoms ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                pStereo[0]->t_parity[iAtom] = atomParity;
                pStereo[0]->nNumber[iAtom] = nAtom1;
                if ( iAtom && !(pStereo[0]->nNumber[iAtom-1] < nAtom1) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        } else {
            for ( p = pStart, iAtom = 0; p < pEnd; iAtom ++, p += (*p==',') ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                if ( !(r = strchr( parity_type, *q) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q+1;
                atomParity = (r - parity_type) + 1;
                pStereo[0]->t_parity[iAtom] = atomParity;
                pStereo[0]->nNumber[iAtom] = nAtom1;
                if ( iAtom && !(pStereo[0]->nNumber[iAtom-1] < nAtom1) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        }
        pStereo[0]->nNumberOfStereoCenters = iAtom;
        /*if ( iAtom ) {*/
            pStereo[0]->nCompInv2Abs = NO_VALUE_INT; /* unknown yet */
        /*}*/

        if ( p != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }

        /* multiplier */
        for ( i = 1; i < mpy_component; i ++ ) {
            ret = CopySegment( pInChI+iComponent+i, pInChI+iComponent, nCpyType, bIso, bIso );
            if ( !ret ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
            ret = CopySegment( pInChI+iComponent+i, pInChI+iComponent, CPY_SP3_M, bIso, bIso );
            if ( !ret ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
            }
            if ( ret < 0 ) {
                goto exit_function;
            }
        }

end_main_cycle:
        iComponent += mpy_component;
        if ( *pEnd ) {
            pStart = pEnd+1;
            continue;
        } else {
            break;
        }
    }
    if ( nNumComponents != iComponent ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentSp2( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[], int state, int *pbAbc )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int i, mpy_component, val;
    int nNumComponents, iComponent, len, iBond;
    AT_NUMB nAtom1, nAtom2;
    int     bondParity;
    const char *p, *q, *t, *pStart, *pEnd, *r;
    int  ret=0;
    INChI *pInChI = pInpInChI[bMobileH];
    INChI *pInChIFrom=NULL;
    /*
    INChI_Stereo *Stereo = NULL;
    INChI_Stereo *StereoOther = NULL;
    */
    INChI_Stereo **pStereo = NULL;

    const char   mult_type[]   = "mnMNe";
    const char   parity_type[] = "-+u?";
    int   bIsoTo, bIsoFrom, nCpyType = CPY_SP2;
    int   bIso = (state==IST_MOBILE_H_ISO_SP2 || state==IST_FIXED_H_ISO_SP2);
    int   base = 10;

    if ( !bIso && state != IST_MOBILE_H_SP2 && state != IST_FIXED_H_SP2 ) {
        return RI_ERR_PROGR; /* program error */
    }
    
    if ( str[0] != 'b' )
        return 0;

    pStart = str+1;
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    if ( !*pStart ) {
        /* creaste empty sp3 segment which means no sp3  */
        for ( iComponent = 0; iComponent < nNumComponents; iComponent ++ ) {
            INChI *pIsoInChI = &pInChI[iComponent];
            pStereo = bIso? &pIsoInChI->StereoIsotopic : &pIsoInChI->Stereo;
            if ( *pStereo && (pStereo[0]->b_parity || pStereo[0]->nNumberOfStereoBonds ||
                 pStereo[0]->nBondAtom1 || pStereo[0]->nBondAtom2 ) ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            /* allocate empty sp3 stereo */
            ret = CopySegment( pIsoInChI, NULL, CPY_SP2, bIso, -1);
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
        ret = nNumComponents+1;
        goto exit_function;
    }

    while( 1 ) {
        
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }

        if ( ((isdigit(*pStart) &&
             0 < (val = (int)inchi_strtol( pStart, &q, 10))) ||
             (q = pStart, val=1))&&
             (t=strchr(mult_type, *q)) && q+1 == pEnd ) {
            /* process the abbreviation */
            ret = 0;
#if (FIX_DALKE_BUGS == 1)
            if ( iComponent + val > nNumComponents ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#endif
            switch( bMobileH ) {
            case TAUT_YES:
                switch( state ) {
                case IST_MOBILE_H_ISO_SP2:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'e' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = -1; /* empty */
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                default:
                    ret = RI_ERR_SYNTAX; 
                    break;
                }
                break;
            case TAUT_NON:
                switch( state ) {
                case IST_FIXED_H_SP2:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to fixed H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 0;
                        bIsoFrom   = 0;
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                case IST_FIXED_H_ISO_SP2:
                    if ( *q == 'm' ) {
                        /* copy from mobile H to fixed isotopic H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'M' ) {
                        /* copy from isotopic mobile H to fixed isotopic H */
                        pInChIFrom = pInpInChI[ALT_TAUT(bMobileH)];
                        bIsoTo     = 1;
                        bIsoFrom   = 1;
                    } else
                    if ( *q == 'n' ) {
                        /* copy from fixed H to fixed isotopic H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = 0;
                    } else
                    if ( *q == 'e' ) {
                        /* copy from mobile H to isotopic mobile H */
                        pInChIFrom = pInChI;
                        bIsoTo     = 1;
                        bIsoFrom   = -1; /* empty */
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                    }
                    break;
                default:
                    ret = RI_ERR_SYNTAX; 
                    break;
                }
                break;

            default:
                ret = RI_ERR_SYNTAX; 
                break;

            }
            if ( ret < 0 ) {
                goto exit_function;
            }
            /* copy */
            for ( i = 0; i < val; i ++ ) {
                ret = CopySegment( pInChI+iComponent+i, pInChIFrom+iComponent+i, nCpyType, bIsoTo, bIsoFrom );
                if ( !ret ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                }
                if ( ret < 0 ) {
                    goto exit_function;
                }
            }
            mpy_component = val;
            goto end_main_cycle;

        } else
        /* regular multiplier */
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
            if ( p != q ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }

            p ++; /* move to the 1st character of the component */
        } else {
            mpy_component = 1;
            p = pStart;
        }
#if (FIX_DALKE_BUGS == 1)
        if ( iComponent + mpy_component > nNumComponents ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
#endif
        pStart = p;
        if ( p < pEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = isupper( UCINT *p)? 1 : 0;
        }
        base = (*pbAbc==1)? ALPHA_BASE : 10;
        if ( *pbAbc == 1 ) {
            /* process the componnt: at1-at2p,at1-at2p,... */
            /* pass 1: find number of stereobonds */
            for ( p = pStart, iBond = 0; p < pEnd; iBond ++ ) {
                /* atoms 1, 2, and parity */
                if ( (nAtom1 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (nAtom2 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (bondParity = (int)inchi_strtol( p, &p, 10),
                      AB_MIN_KNOWN_PARITY <= bondParity && bondParity <= AB_MAX_KNOWN_PARITY) ) {
                    ; /* okay */
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( nAtom1 <= nAtom2   ||
                     nAtom1 > pInChI[iComponent].nNumberOfAtoms ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        } else {
            /* process the componnt: at1-at2p,at1-at2p,... */
            /* pass 1: find number of stereobonds */
            for ( p = pStart, iBond = 0; p < pEnd; iBond ++, p += (*p==',') ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                if ( *q != '-' ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q+1;
                nAtom2 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                if ( !nAtom1 || !nAtom2 ||
                     nAtom1 <= nAtom2   ||
                     nAtom1 > pInChI[iComponent].nNumberOfAtoms ||
                     !(r = strchr( parity_type, *q) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q+1;
            }
        }

        if ( p != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        len = iBond;

        /* memory allocation */

        pStereo = bIso? &pInChI[iComponent].StereoIsotopic : &pInChI[iComponent].Stereo;

        if ( !*pStereo ) {
            if ( !(*pStereo = (INChI_Stereo *) inchi_calloc( 1, sizeof(**pStereo) ) ) ) {
                ret = RI_ERR_ALLOC; /* memory allocation failed */
                goto exit_function;
            }
        }
        if ( pStereo[0]->b_parity || pStereo[0]->nNumberOfStereoBonds ||
             pStereo[0]->nBondAtom1 || pStereo[0]->nBondAtom2 ) {
            ret = RI_ERR_SYNTAX; /* syntax error: bonds have already been allocated */
            goto exit_function;
        }
        /* allocate sp2 stereo */
        if ( !(pStereo[0]->b_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pStereo[0]->b_parity[0]) ) ) ||
             !(pStereo[0]->nBondAtom1 = (AT_NUMB *)inchi_calloc( len+1, sizeof(pStereo[0]->nBondAtom1[0]) ) ) ||
             !(pStereo[0]->nBondAtom2 = (AT_NUMB *)inchi_calloc( len+1, sizeof(pStereo[0]->nBondAtom2[0]) ) ) ) {
            /* cleanup */
            if ( pStereo[0]->b_parity ) {
                INCHI_HEAPCHK
                inchi_free( pStereo[0]->b_parity );
                pStereo[0]->b_parity = NULL;
            }
            if ( pStereo[0]->nBondAtom1 ) {
                INCHI_HEAPCHK
                inchi_free( pStereo[0]->nBondAtom1 );
                pStereo[0]->nBondAtom1 = NULL;
            }
            if ( pStereo[0]->nBondAtom2 ) {
                INCHI_HEAPCHK
                inchi_free( pStereo[0]->nBondAtom2 );
                pStereo[0]->nBondAtom2 = NULL;
            }
            INCHI_HEAPCHK
            ret = RI_ERR_ALLOC; /* memory allocation failed */
            goto exit_function;
        }

        /* pass 2: store stereobonds */
        if ( *pbAbc == 1 ) {
            for ( p = pStart, iBond = 0; p < pEnd; iBond ++ ) {
                if ( (nAtom1 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (nAtom2 = (AT_NUMB)inchi_strtol( p, &p, base ) ) &&
                      (bondParity = (int)inchi_strtol( p, &p, 10),
                      AB_MIN_KNOWN_PARITY <= bondParity && bondParity <= AB_MAX_KNOWN_PARITY) ) {
                    ; /* okay */
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                pStereo[0]->b_parity[iBond] = bondParity;
                pStereo[0]->nBondAtom1[iBond] = nAtom1;
                pStereo[0]->nBondAtom2[iBond] = nAtom2;

                if ( iBond && 
                     !(pStereo[0]->nBondAtom1[iBond-1] <  nAtom1 ||
                       (pStereo[0]->nBondAtom1[iBond-1] == nAtom1 &&
                       pStereo[0]->nBondAtom2[iBond-1] <  nAtom2) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error: wrong bond order */
                    goto exit_function;
                }
            }
        } else {
            for ( p = pStart, iBond = 0; p < pEnd; iBond ++, p += (*p==',') ) {
                nAtom1 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                if ( *q != '-' ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q+1;
                nAtom2 = (AT_NUMB)inchi_strtol( p, &q, 10 );
                if ( !(r = strchr( parity_type, *q) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q+1;
                bondParity = (r - parity_type) + 1;
                pStereo[0]->b_parity[iBond] = bondParity;
                pStereo[0]->nBondAtom1[iBond] = nAtom1;
                pStereo[0]->nBondAtom2[iBond] = nAtom2;

                if ( iBond && 
                     !(pStereo[0]->nBondAtom1[iBond-1] <  nAtom1 ||
                       (pStereo[0]->nBondAtom1[iBond-1] == nAtom1 &&
                       pStereo[0]->nBondAtom2[iBond-1] <  nAtom2) ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error: wrong bond order */
                    goto exit_function;
                }
            }
        }
        pStereo[0]->nNumberOfStereoBonds = iBond;

        if ( p != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }

        /* multiplier */
        for ( i = 1; i < mpy_component; i ++ ) {
            ret = CopySegment( pInChI+iComponent+i, pInChI+iComponent, nCpyType, bIso, bIso );
            if ( ret < 0 ) {
                goto exit_function;
            }
        }

end_main_cycle:
        iComponent += mpy_component;
        if ( *pEnd ) {
            pStart = pEnd+1;
            continue;
        } else {
            break;
        }
    }
    if ( nNumComponents != iComponent ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentProtons( const char *str, int bMobileH, REM_PROTONS nNumProtons[], int ppnNumComponents[] )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int val;
    const char *q, *pStart, *pEnd;
    int  ret;
    
    if ( str[0] != 'p' )
        return 0;

    pStart = str+1;

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }

        if ( pStart[0] == '+' && isdigit( UCINT pStart[1] ) ) {
            val = (int)inchi_strtol( pStart+1, &q, 10 );
        } else
        if ( pStart[0] == '-' && isdigit( UCINT pStart[1] ) ) {
            val = -(int)inchi_strtol( pStart+1, &q, 10 );
        } else {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        if ( !val ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        nNumProtons[bMobileH].nNumRemovedProtons = val;
        if ( *pEnd || q != pEnd ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        } else {
            break;
        }
    }
    ret = 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentCharge( const char *str, int bMobileH, INChI *pInpInChI[], int ppnNumComponents[] )
{
    /* Pass 1: count bonds and find actual numbers of  atom */
    int i, mpy_component, val;
    int nNumComponents, iComponent;
    const char *p, *q, *t, *pStart, *pEnd;
    int  ret;
    INChI *pInChI = pInpInChI[bMobileH];
    const char   mult_type[] = "mnMNe";
    
    if ( str[0] != 'q' ) {
        return 0;
    }

    pStart = str+1;
    iComponent = 0;
    nNumComponents = ppnNumComponents[bMobileH];

    if ( !*pStart && bMobileH == TAUT_NON ) {
        for ( i = 0; i < nNumComponents; i ++ ) {
            pInChI[i].nTotalCharge = NO_VALUE_INT;
        }
        return nNumComponents+1;
    }

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }

        if ( ((isdigit(UCINT *pStart) &&
             0 < (val = (int)inchi_strtol( pStart, &q, 10))) ||
             (q = pStart, val=1) )&&
             (t=strchr(mult_type, *q)) && q+1 == pEnd ) {
            /* process the abbreviation */

            switch( bMobileH ) {

            case TAUT_YES:
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;

            case TAUT_NON:
                if ( *q != 'm' ||
                    iComponent + val > nNumComponents ||
                    iComponent + val > ppnNumComponents[TAUT_YES] ) {

                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                for ( i = 0; i < val; i ++ ) {
                    /* avoid 0 which means "omitted" */
                    pInChI[iComponent+i].nTotalCharge = pInpInChI[TAUT_YES][iComponent+i].nTotalCharge?
                                                        pInpInChI[TAUT_YES][iComponent+i].nTotalCharge :
                                                        NO_VALUE_INT;
                }
                mpy_component = val;
                goto end_main_cycle;

            default:
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
        }else
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
            if ( p != q ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            p ++;
        } else {
            mpy_component = 1;
            p = pStart;
        }
#if ( FIX_DALKE_BUGS == 1 )
        if ( mpy_component + iComponent > nNumComponents || mpy_component <= 0 ) {
            ret = RI_ERR_SYNTAX; /* syntax error: too many components in charge layer */
            goto exit_function;
        }
#endif
        pStart = p;

        if ( pStart < pEnd ) {
            if ( pStart[0] == '+' && isdigit( UCINT pStart[1] ) ) {
                val = (int)inchi_strtol( pStart+1, &q, 10 );
                pStart = q;
            } else
            if ( pStart[0] == '-' && isdigit( UCINT pStart[1] ) ) {
                val = -(int)inchi_strtol( pStart+1, &q, 10 );
                pStart = q;
            } else {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#if ( FIX_DALKE_BUGS == 1 )
            if ( val < -256 || val > 256 ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
#endif
            if ( !val ) {
                if ( pStart != pEnd ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( bMobileH == TAUT_NON ) {
                    val = NO_VALUE_INT;  /* avoid 0 which means "omitted" */
                }
            }
        } else {
            val = NO_VALUE_INT;
        }
        for ( i = 0; i < mpy_component; i ++ ) {
            pInChI[iComponent+i].nTotalCharge = val;
        }

end_main_cycle:
        iComponent += mpy_component;
        if ( *pEnd ) {
            pStart = pEnd+1;
            continue;
        } else {
            break;
        }
    }

    if ( nNumComponents != iComponent ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    return ret;   

}
/****************************************************************************************/
int ParseSegmentMobileH( const char *str, int bMobileH, INChI *pInpInChI[], int pnNumComponents[], int *pbAbc )
{
#define nNum_H( ICOMPONENT ) ((bMobileH==TAUT_YES)? pInChI[ICOMPONENT].nNum_H : pInChI[ICOMPONENT].nNum_H_fixed)    
    /* Pass 1: count bonds and find actual numbers of  atom */
    int i, mpy_component, num_H, num_Minus, val, num_Atoms, numCtAtoms, tg_alloc_len, len, len2;
    int num_H_component, num_H_formula, num_taut_H_component, num_H_InChI, ret2;
    int nNumComponents, iComponent, nNumBonds, lenTautomer, tg_pos_Tautomer, iTGroup;
    const char *p, *q, *h, *t, *p1, *pTaut, *pStart, *pEnd;
    AT_NUMB curAtom, nxtAtom;
    int  num_open, state, ret, nAltMobileH = ALT_TAUT(bMobileH);
    INChI *pInChI = pInpInChI[bMobileH];
    INChI *pAltInChI = pInpInChI[nAltMobileH];
    int  base = 10;
    
    num_H = -999;          /* impossible value */
    num_Minus = -999;      /* impossible value */
    tg_pos_Tautomer = -999; /* impossible value */

    /* number of immobile H is always allocated; immobile H are present in M layer only */
    nNumComponents = pnNumComponents[bMobileH];
    for ( i = 0; i < nNumComponents; i ++ ) {
        len = pInChI[i].nNumberOfAtoms;
        if ( bMobileH == TAUT_NON && i < pnNumComponents[nAltMobileH] ) {
            if ( len < pAltInChI[i].nNumberOfAtoms ) {
                len = pAltInChI[i].nNumberOfAtoms;
                if ( pInChI[i].nNum_H ) {
                    inchi_free( pInChI[i].nNum_H );
                    pInChI[i].nNum_H = NULL;
                }
            }
        }
        len ++;
        if ( !pInChI[i].nNum_H && /* allocate immobile H segment if it has not been allocated yet */ 
             !(pInChI[i].nNum_H = (S_CHAR *)inchi_calloc( len, sizeof(pInChI[0].nNum_H[0]) )) ) {
            ret = RI_ERR_ALLOC; /* allocation error */
            goto exit_function;
        }
        /* copy immobile H from Mobile-H layer to Fixed-H layer */
        if ( bMobileH == TAUT_NON && i < pnNumComponents[nAltMobileH] ) {
            memcpy( pInChI[i].nNum_H, pAltInChI[i].nNum_H, (len-1) * sizeof(pInChI[0].nNum_H[0]) );
        }
    }

    if ( str[0] != 'h' )
        return 0;

    /* Read Hydrogen info in 1 pass */
    pStart = str+1;
    iComponent = 0;
    nNumComponents = pnNumComponents[bMobileH];

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
#if ( FIX_DALKE_BUGS == 1 )
            if ( p != q || !isdigit(UCINT *pStart) ) /* prevent non-positive multipliers */
#else
            if ( p != q )
#endif
            {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            p ++;
        } else {
            mpy_component = 1;
            p = pStart;
        }
        pStart = p;
        /* Pass 1.1 parse a component */
        num_open = 0;
        state='\0';   /* initial state */
        nNumBonds = 0;
        curAtom   = 0;
        numCtAtoms = pInChI[iComponent].nNumberOfAtoms;
        if ( bMobileH == TAUT_NON && iComponent < pnNumComponents[nAltMobileH] ) {
            numCtAtoms = pAltInChI[iComponent].nNumberOfAtoms;
        }
        
        if ( p < pEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = (*p == ',' || isupper( UCINT *p))? 1 : 0;
        }
        base = (*pbAbc==1)? ALPHA_BASE : 10;

        /* immobile H */
        t = pTaut = (*pbAbc==1)? strchr( p, ',' ) : strchr( p, '(' ); /* locate the first tautomer group character */
        
        if ( t && bMobileH == TAUT_NON ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
        if ( !pTaut || pTaut > pEnd ) {
            pTaut = pEnd;
            t = NULL; /* found no tautomeric group for this component */
        }
        for ( i = 0; i < mpy_component; i ++ ) {
            if ( bMobileH == TAUT_NON ) {
                /* allocate nNum_H_fixed */
                if ( pInChI[iComponent+i].nNum_H_fixed ) {
                    ret = RI_ERR_PROGR; /* program error */
                    goto exit_function;
                }
                if ( iComponent+i < pnNumComponents[nAltMobileH] ) {
                    len = inchi_max(pInChI[iComponent+i].nNumberOfAtoms, pAltInChI[iComponent+i].nNumberOfAtoms)+1;
                } else {
                    len = pInChI[iComponent+i].nNumberOfAtoms + 1;
                }
                pInChI[iComponent+i].nNum_H_fixed = (S_CHAR *)inchi_calloc( len, sizeof(pInChI[0].nNum_H_fixed[0]) );
                if ( !pInChI[iComponent+i].nNum_H_fixed ) {
                    ret = RI_ERR_ALLOC; /* allocation error */
                    goto exit_function;
                }
                /* compare nAtom */
                if ( iComponent+i < pnNumComponents[nAltMobileH] ) {
                    len2 = inchi_min(pInChI[iComponent+i].nNumberOfAtoms, pAltInChI[iComponent+i].nNumberOfAtoms);
                    if ( pInChI[iComponent+i].nAtom && len2 ) {
                        /* check */
                        if ( memcmp( pInChI[iComponent+i].nAtom, pAltInChI[iComponent+i].nAtom, len2 * sizeof(pInChI[0].nAtom[0]) ) ) {
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    }
                    /* allocate and copy atom if bridging H are present */
                    if ( pInChI[iComponent+i].nNumberOfAtoms < pAltInChI[iComponent+i].nNumberOfAtoms ) {
                        if ( pInChI[iComponent+i].nAtom )
                            inchi_free( pInChI[iComponent+i].nAtom );
                        if ( !(pInChI[iComponent+i].nAtom = (U_CHAR *)inchi_calloc( len, sizeof(pInChI[0].nAtom[0]) ) ) ) {
                            ret = RI_ERR_ALLOC; /* allocation error */
                            goto exit_function;
                        }
                        if ( len > 1 ) {
                            memcpy( pInChI[iComponent+i].nAtom, pAltInChI[iComponent+i].nAtom, (len-1) * sizeof(pInChI[0].nAtom[0]) );
                        }
                        /* correct number of atoms including bridging H */
                        pInChI[iComponent+i].nNumberOfAtoms = pAltInChI[iComponent+i].nNumberOfAtoms;
                    }
                }
            }
        }
        
        if ( *pbAbc == 1 ) {
            /* read numbers of H: XnYn... or XYn... */
            p = pStart;
            tg_alloc_len = 0;
            num_H_component = num_taut_H_component = 0;
            while ( p < pTaut ) {
                /* syntax check: atom number */
                if ( !*p || !isupper( UCINT *p ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( (curAtom = nxtAtom = (int)inchi_strtol( p, &q, base )) ) {
                    p = q;
                    if ( isupper( UCINT *p ) ) {
                        nxtAtom = (int)inchi_strtol( p, &q, base );
                        p = q;
                    }
                }
                if ( curAtom > nxtAtom || nxtAtom > numCtAtoms || p > pTaut ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                /* number of H, may be negative */
                if ( !(num_H = (int)inchi_strtol( p, &q, 10 )) || q > pTaut ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                p = q;
                /* set number of H */
                for ( i = curAtom; i <= nxtAtom; i ++ ) {
                    nNum_H(iComponent)[i-1] = num_H;
                    num_H_component += num_H;
                }
            }
            if ( p != pTaut ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
        } else {
            /* read numbers of H: 1-2,3H2,4,5H3 */
            p = pStart;
            tg_alloc_len = 0;
            num_H_component = num_taut_H_component = 0;
            while ( p < pTaut ) {
                /* syntax check: atom number */
                if ( !*p || !isdigit( UCINT *p ) ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                /* number of H */
                h = p + strcspn( p, "Hh" );
                /*h = strchr( p, 'H' );*/
                if ( !*h || h >= pTaut ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                
                    /*
                    p = pTaut;
                    h = NULL;
                    break; */ /* no more H found */
                }
                num_H = (*h == 'H')? 1 : (*h == 'h')? -1 : 0;
                if ( !num_H ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( h[1] && isdigit( UCINT h[1] ) ) {
                    num_H *= (int)inchi_strtol( h+1, &p1, 10 );
                } else {
                    p1 = h+1; /* next set of immobile H */
                }
                if ( *p1 == ',' ) {
                    p1 ++;  /* next H-subsegment; otherwise (H or ; or end of the segment */
                }
                /* list of atoms that have num_H */
                while ( p < h ) {
                    if ( !*p || !isdigit( UCINT *p ) ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    nxtAtom = curAtom = (int)inchi_strtol( p, &q, 10 );
                    if ( *q == '-' ) {
                        nxtAtom = (int)inchi_strtol( q+1, &q, 10 );
                    }
                    /* consitency check */
                    if ( !curAtom || curAtom > numCtAtoms ||
                         nxtAtom < curAtom || nxtAtom > numCtAtoms ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    /* set number of H */
                    for ( i = curAtom; i <= nxtAtom; i ++ ) {
                        nNum_H(iComponent)[i-1] = num_H;
                        num_H_component += num_H;
                    }
                    /* move to the next atom number if any */
                    p = q;
                    if ( *p == ',' ) {
                        p ++;
                    }
                }

                if ( p == h ) {
                    p = p1;
                } else
                if ( p == pTaut ) {
                    break;
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        }

        INCHI_HEAPCHK
        /* ) -> (, H, N, [-, N,], AtNum,... AtNum) */ 
        lenTautomer = 0;
        if ( (p = t) ) {
            if ( *pbAbc == 1 ) {
                /* tautomeric groups: pass 1 */
                iTGroup = 0;
                state   = ')';  /* init as if the prev. t-group just ended */
                num_Atoms = 0;
                /* Tautomeric info storage */
                /* NumGroups; ((NumAt+2, NumH, Num(-), At1..AtNumAt),...); {INCHI_T_NUM_MOVABLE = 2} */
                /* Allocated length: [5*nNumberOfAtoms/2+1], see Alloc_INChI(...) */
                if ( *p == ',' ) { /* start t-group */
                    p ++;
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                while ( p < pEnd ) {
                    /* start t-group */
                    if ( !isdigit( UCINT *p ) || !(num_H = (int)inchi_strtol( p, &q, 10 ) ) || q > pEnd ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    p = q;
                    num_Minus = 0;
                    if ( *p == '-' ) {
                        p ++;
                        if ( isdigit( UCINT *p ) ) {
                            num_Minus = (int)inchi_strtol( p, &q, 10 );
                            p = q;
                        } else {
                            num_Minus = 1;
                        }
                    }
                    if ( p >= pEnd ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    if ( !tg_alloc_len ) {
                        /* 
                           --- header ---
                           [num_t_groups]
                           --- one t-group: ---
                           [len=group length no including this value]
                           [num_H]
                           [num_(-)]
                           [Endpoint(1),...,Endpoint(len-2)]
                           --- next t-group ---
                           ...

                           Max. size = 1 + 3*max_num_t_groups + max_num_endpoints

                           max_num_t_groups  = num_at/2
                           max_num_endpoints = num_at

                           Max. size = 1 + 3*(num_at/2) + num_at = 1 + (5*num_at)/2
                           5 = 3 + INCHI_T_NUM_MOVABLE = 3 + num_types_of_attachments

                           This does not include zero termination!
                       
                        */
                        tg_alloc_len = ((3+INCHI_T_NUM_MOVABLE)*pInChI[iComponent].nNumberOfAtoms)/2+1;
                        for ( i = 0; i < mpy_component; i ++ ) {
                            pInChI[iComponent+i].nTautomer = (AT_NUMB*)inchi_calloc( tg_alloc_len+1, sizeof(pInChI->nTautomer[0]));
                            if ( !pInChI[iComponent+i].nTautomer ) {
                                ret = RI_ERR_ALLOC; /* allocation error */
                                goto exit_function;
                            }
                            pInChI[iComponent+i].lenTautomer = 0;
                        }
                        tg_pos_Tautomer = 1; /* number atoms (NumAt+2) position */
                    } else {
                        /* next t-group */
                        tg_pos_Tautomer = lenTautomer;
                    }
                    if ( tg_pos_Tautomer+3 >= tg_alloc_len ) {
                        ret = RI_ERR_PROGR; /* wrong tautomer array length */
                        goto exit_function;
                    }
                    pInChI[iComponent].nTautomer[tg_pos_Tautomer+1] = num_H;
                    pInChI[iComponent].nTautomer[tg_pos_Tautomer+2] = num_Minus;
                    lenTautomer                = tg_pos_Tautomer+3; /* first atom number position */
                    num_taut_H_component += num_H;

                    while ( p < pEnd && isupper( UCINT *p) ) {
                        /* read list of tautomeric atoms */
                        val = (int)inchi_strtol( p, &q, base );
                        if ( lenTautomer >= tg_alloc_len || val > numCtAtoms ) {
                            ret = RI_ERR_PROGR; /* wrong tautomer array length */
                            goto exit_function;
                        }
                        num_Atoms ++;
                        pInChI[iComponent].nTautomer[lenTautomer ++] = val;
                        p = q;
                    }
                    if ( !num_Atoms || (p < pEnd && !isdigit( UCINT *p))  ) {
                        ret = RI_ERR_PROGR; /* wrong tautomer array length */
                        goto exit_function;
                    }
                    iTGroup ++;
                    pInChI[iComponent].nTautomer[tg_pos_Tautomer] = lenTautomer - tg_pos_Tautomer - 1; /* length of the rest of the t-group */
                    pInChI[iComponent].lenTautomer = lenTautomer;
                }
                if ( !iTGroup || p != pEnd ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                pInChI[iComponent].nTautomer[0] = iTGroup;

            } else {
                /* tautomeric groups: pass 1 */
                iTGroup = 0;
                state   = ')';  /* init as if the prev. t-group just ended */
                num_Atoms = 0;
                /* Tautomeric info storage */
                /* NumGroups; ((NumAt+2, NumH, Num(-), At1..AtNumAt),...); {INCHI_T_NUM_MOVABLE = 2} */
                /* Allocated length: [5*nNumberOfAtoms/2+1], see Alloc_INChI(...) */
                while ( p < pEnd ) {
                    /* t-group */
                    switch ( *p ) {
                    case '(': /* start t-group */
                        switch ( state ) {
                        case ')':
                            state = *p ++;
                            num_H = 0;
                            num_Minus = 0;
                            continue;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    case ')': /* end t-group */
                        switch ( state ) {
                        case 'A': /* previuos was atom number */
                            if ( !tg_alloc_len ) {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                            iTGroup ++;
                            state = *p ++;
                            pInChI[iComponent].nTautomer[tg_pos_Tautomer] = lenTautomer - tg_pos_Tautomer - 1; /* length of the rest of the t-group */
                            pInChI[iComponent].lenTautomer = lenTautomer;
                            continue;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    case 'H': /* number of H */
                        switch ( state ) {
                        case '(':
                            state = *p ++;
                            num_H = 1;
                            continue;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    case '-':  /* number of (-) */
                        switch ( state ) {
                        case 'N': /* previous was number of H */
                        case 'H': /* previous was H */
                            state = *p ++;
                            num_Minus = 1;
                            continue;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    case ',':
                        switch ( state ) {
                        case 'N': /* previous was number of H */
                        case 'H': /* previous was H */
                        case '-': /* previuos was - */
                        case 'M': /* previous was number of (-) */
                            /* the next must be the first tautomeric atom number; save num_H & num_Minus */
                            if ( num_H <= 0 && num_Minus <= 0 ) {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                            if ( !tg_alloc_len ) {
                                /* 
                                   --- header ---
                                   [num_t_groups]
                                   --- one t-group: ---
                                   [len=group length no including this value]
                                   [num_H]
                                   [num_(-)]
                                   [Endpoint(1),...,Endpoint(len-2)]
                                   --- next t-group ---
                                   ...

                                   Max. size = 1 + 3*max_num_t_groups + max_num_endpoints

                                   max_num_t_groups  = num_at/2
                                   max_num_endpoints = num_at

                                   Max. size = 1 + 3*(num_at/2) + num_at = 1 + (5*num_at)/2
                                   5 = 3 + INCHI_T_NUM_MOVABLE = 3 + num_types_of_attachments

                                   This does not include zero termination!
                               
                                */
                                tg_alloc_len = ((3+INCHI_T_NUM_MOVABLE)*pInChI[iComponent].nNumberOfAtoms)/2+1;
                                for ( i = 0; i < mpy_component; i ++ ) {
                                    pInChI[iComponent+i].nTautomer = (AT_NUMB*)inchi_calloc( tg_alloc_len+1, sizeof(pInChI->nTautomer[0]));
                                    if ( !pInChI[iComponent+i].nTautomer ) {
                                        ret = RI_ERR_ALLOC; /* allocation error */
                                        goto exit_function;
                                    }
                                    pInChI[iComponent+i].lenTautomer = 0;
                                }
                                tg_pos_Tautomer = 1; /* number atoms (NumAt+2) position */
                            } else {
                                /* next t-group */
                                tg_pos_Tautomer = lenTautomer;
                            }
                            if ( tg_pos_Tautomer+3 >= tg_alloc_len ) {
                                ret = RI_ERR_PROGR; /* wrong tautomer array length */
                                goto exit_function;
                            }
                            pInChI[iComponent].nTautomer[tg_pos_Tautomer+1] = num_H;
                            pInChI[iComponent].nTautomer[tg_pos_Tautomer+2] = num_Minus;
                            lenTautomer                = tg_pos_Tautomer+3; /* first atom number position */
                            num_taut_H_component += num_H;
                            state = *p ++;
                            continue;
                        case 'A':  /* previuos was atom number */
                            state = *p ++;
                            continue;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    default:
                        if ( isdigit( UCINT *p ) ) {
                            val = (int)inchi_strtol( p, &q, 10 );
                            if ( val <= 0 ) {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                            p = q;
                            switch( state ) {
                            case 'H':
                                num_H = val;
                                state = 'N';
                                continue;
                            case '-':
                                num_Minus = val;
                                state = 'M';
                                continue;
                            case ',':
                                if ( lenTautomer >= tg_alloc_len || val > numCtAtoms ) {
                                    ret = RI_ERR_PROGR; /* wrong tautomer array length */
                                    goto exit_function;
                                }
                                num_Atoms ++;
                                pInChI[iComponent].nTautomer[lenTautomer ++] = val;
                                state = 'A';
                                continue;
                            default:
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;

                            }
                        }
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                }
                if ( !iTGroup || state != ')' ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                pInChI[iComponent].nTautomer[0] = iTGroup;
            }
        }
        /* check num_H in components; for bMobileH=TAUT_NON, pInChI->nNum_H_fixed[] has not been added to pInChI->nNum_H[] yet */
        if ( 0 > ( ret2 = GetInChIFormulaNumH( pInChI+iComponent, &num_H_formula) ) ||
             0 > ( ret2 = GetInChINumH( pInChI+iComponent, &num_H_InChI ) ) ) {
            ret = ret2;
            goto exit_function;
        }
        if ( num_H_formula != num_H_InChI + (bMobileH==TAUT_NON? num_H_component:0) ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }

        /* duplicate according to multipolication */
        for ( i = 1; i < mpy_component; i ++ ) {
            memcpy( nNum_H(iComponent+i), nNum_H(iComponent), pInChI[iComponent+i].nNumberOfAtoms * sizeof(nNum_H(0)[0]) );
            /*
            memcpy( pInChI[iComponent+i].nNum_H, pInChI[iComponent].nNum_H,
                    pInChI[iComponent+i].nNumberOfAtoms * sizeof(pInChI[0].nNum_H[0]) );
            */
            if ( pInChI[iComponent+i].nTautomer && pInChI[iComponent].nTautomer && pInChI[iComponent].lenTautomer ) {
                memcpy( pInChI[iComponent+i].nTautomer, pInChI[iComponent].nTautomer,
                        pInChI[iComponent].lenTautomer * sizeof(pInChI[0].nTautomer[0]) );
                pInChI[iComponent+i].lenTautomer = pInChI[iComponent].lenTautomer;
            }
            /* check num_H in components */
            if ( 0 > ( ret2 = GetInChIFormulaNumH( pInChI+iComponent+i, &num_H_formula) ) ||
                 0 > ( ret2 = GetInChINumH( pInChI+iComponent+i, &num_H_InChI ) ) ) {
                ret = ret2;
                goto exit_function;
            }
            if ( num_H_formula != num_H_InChI + (bMobileH==TAUT_NON? num_H_component:0) ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
        }

        /* prepare for the next component */
        iComponent += i;
        if ( *pEnd ) {
#if (FIX_DALKE_BUGS == 1)
            /* prevent crash on extra trailing ';' */
            if ( iComponent >= nNumComponents ) {
                ret = RI_ERR_SYNTAX; /* syntax error: extra component */
                goto exit_function;
            }
#endif
            pStart = pEnd+1;
        } else
            break;
    }
    if ( nNumComponents != iComponent ) {
        ret = RI_ERR_SYNTAX; /* syntax error */
        goto exit_function;
    }
    ret = iComponent + 1;

exit_function:
    INCHI_HEAPCHK
    return ret;   

}

/****************************************************************************************/
int ParseSegmentConnections( const char *str, int bMobileH, INChI **pInpInChI, int *pnNumComponents, int *pbAbc )
{
#define LAST_AT_LEN 256
    /* Pass 1: count bonds and find actual numbers of  atom */
    int i, j, k, m, c, mpy_component;
    int nNumComponents, iComponent, nNumAtoms, nNumBonds, lenConnTable, iBond;
    const char *p, *q, *pStart, *pEnd;
    AT_NUMB last_atom[LAST_AT_LEN], curAtom, maxAtom;
    int  num_open, state, ret, base;
    INChI *pInChI = *pInpInChI;
    LINKED_BONDS LB;
    LINKED_BONDS *pLB = &LB; /* a list of linked lists of bonds, for each atom */
    AT_NUMB neighbor[MAXVAL];
    int bPrevVersion = -1;
    
    iComponent = 0;
    if ( str[0] != 'c' ) {
        if ( !pInChI && !*pnNumComponents ) {
            int lenFormula = 1;
            /* component has no formula; allocate InChI */
            lenConnTable   = 0;
            nNumComponents = 1;
            /* allocate InChI */
            if ( !(pInChI = *pInpInChI = (INChI *)inchi_calloc( nNumComponents, sizeof(INChI) ) ) ) {
                return RI_ERR_ALLOC; /* alloc failure */
            }
            /* allocate empty formula */
            pInChI[iComponent].szHillFormula = (char *)inchi_calloc( lenFormula+1, sizeof(pInChI[0].szHillFormula[0]) );
            if ( !pInChI[iComponent].szHillFormula ) {
                ret = RI_ERR_ALLOC; /* allocation failure */
                goto exit_function;
            }
            /* allocate empty connection table */
            pInChI[iComponent].nConnTable = (AT_NUMB *)inchi_calloc( lenConnTable+1, sizeof(pInChI[0].nConnTable[0]) );
            if ( !pInChI[iComponent].nConnTable ) {
                ret = RI_ERR_ALLOC; /* allocation failure */
                goto exit_function;
            }
            pInChI[iComponent].lenConnTable = lenConnTable;
            *pnNumComponents = nNumComponents;
        } else {
            lenConnTable   = 1;
            nNumComponents = *pnNumComponents;
            for ( i = 0; i < nNumComponents; i ++ ) {
                /* allocate 1 atom connection table */
                if ( pInChI[i].nConnTable ) {
                    inchi_free( pInChI[i].nConnTable );
                }
                pInChI[i].nConnTable = (AT_NUMB *)inchi_calloc( lenConnTable+1, sizeof(pInChI[0].nConnTable[0]) );
                if ( !pInChI[i].nConnTable ) {
                    ret = RI_ERR_ALLOC; /* allocation failure */
                    goto exit_function;
                }
                pInChI[i].nConnTable[0] = 1;
                pInChI[i].lenConnTable = lenConnTable;
            }
        }
        return 0;
    }

    /* Pass 1. Re-Count atoms, count bonds */
    pStart = str+1;
    nNumComponents = *pnNumComponents;
#if (FIX_DALKE_BUGS == 1)
    /* prevent crash on too many components */
    if ( nNumComponents > MAX_ATOMS ) {
        ret = RI_ERR_SYNTAX; /* syntax error: extra component */
        goto exit_function;
    }
#endif
    memset( pLB, 0, sizeof(pLB[0]) );

    while( 1 ) {
        /* cycle over components */
        if ( !(pEnd = strchr( pStart, ';' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        if ( (p = strchr(pStart, '*')) && p < pEnd ) {
            mpy_component = (int)inchi_strtol( pStart, &q, 10 );
            if ( p != q 
#if (FIX_DALKE_BUGS == 1)
                || !isdigit( UCINT *pStart )
#endif
                ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            p ++;
        } else {
            mpy_component = 1;
            p = pStart;
        }
#if (FIX_DALKE_BUGS == 1)
        if ( iComponent + mpy_component > MAX_ATOMS ) {
            ret = RI_ERR_SYNTAX; /* syntax error */
            goto exit_function;
        }
#endif
        pStart = p;
        /* Pass 1.1 parse a component */
        num_open = 0;
        memset( last_atom, 0, sizeof(last_atom) );
        state='\0';   /* initial state */
        maxAtom = 0;
        nNumBonds = 0;
        curAtom   = 0;
        if ( p < pEnd && *pbAbc == -1 ) {
            /* check if compressed InChI */
            *pbAbc = isupper( UCINT *p)? 1 : 0;
        }
        base = *pbAbc? ALPHA_BASE : 10;

        if ( *pbAbc == 1 ) {
            nNumAtoms = 1;
            while ( p < pEnd ) {
                if ( *p == '-' ) {
                    if ( bPrevVersion == -1 ) {
                        /* previous InChI version */
                        bPrevVersion = 1;
                    } else
                    if ( bPrevVersion != 1 ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    nNumAtoms --;
                    p ++;
                }
                if ( isdigit( UCINT *p ) ) {
                    if ( bPrevVersion == -1 ) {
                        /* curreny InChI, version 1 */
                        bPrevVersion = 0;
                    } else
                    if ( bPrevVersion != 0 ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    nNumAtoms -= inchi_strtol( p, &p, 10 ); /* bypass digits */
                }
                if ( *p != '-' &&  ( curAtom = (AT_NUMB)inchi_strtol( p, &q, base ) ) ) {
                    nNumAtoms ++;
                    nNumBonds ++;
                    p = q;
                    if ( maxAtom < curAtom )
                        maxAtom = curAtom;
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
            if ( maxAtom < nNumAtoms && nNumBonds ) {
                maxAtom = nNumAtoms;
            }
        } else {
            while ( p < pEnd ) {
                /* atom number */
                c = UCINT *p ++;
                switch ( c ) {
                case '(':
                case ')':
                case ',':
                case '-':
                    if ( state != 'N' ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    state = c;
                    num_open += (c=='(') - (c==')');
                    if ( num_open < 0 ) {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    break;
                default:
                    if ( isdigit( c ) && (curAtom = (AT_NUMB)inchi_strtol( p-1, &q, 10 )) ) {
                        p = q;
                        switch( state ) {
                        case '(':
                        case ')':
                        case ',':
                        case '-':
                            nNumBonds ++;
                        case '\0':
                            if ( maxAtom < curAtom )
                                maxAtom = curAtom;
                            state = 'N';
                            break;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    break;
                }
            }
            if ( num_open ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
                /* syntax error: parentheses do not match */
            }
        }
        /* Save the results and allocate memory */
        nNumAtoms = (int)maxAtom; /* 0 if empty connection table and no bonds present */
        lenConnTable = nNumAtoms + nNumBonds;
        /* connection table format: At1[,Neigh11,Neigh12,...],At2[,Neigh21,Neigh22,...],AtN[NeighN1,NeighN2,...] */
        /* where AtK > NeighK1 > NeighK2,...; At(K) < At(K+1); the length = num.atoms + num.bonds */
        for ( i = 0; i < mpy_component; i ++ ) {
            /* check number of atoms: the difference may be due to bridging H */
            if ( (j = pInChI[iComponent+i].nNumberOfAtoms) < nNumAtoms ) {
                /* reallocate */
                U_CHAR *nAtomTmp = (U_CHAR *) inchi_malloc( nNumAtoms + 1 );
                if ( !nAtomTmp ) {
                    ret = RI_ERR_ALLOC; /* allocation failure */
                    goto exit_function;
                }
                memcpy( nAtomTmp, pInChI[iComponent+i].nAtom, sizeof(nAtomTmp[0])*j);
                while ( j < nNumAtoms ) {
                    nAtomTmp[j ++] = EL_NUMBER_H; /* bridging H */
                }
                nAtomTmp[j] = '\0';
                INCHI_HEAPCHK
                if ( pInChI[iComponent+i].nAtom ) {
                    inchi_free( pInChI[iComponent+i].nAtom );
                }
                pInChI[iComponent+i].nAtom = nAtomTmp;
                pInChI[iComponent+i].nNumberOfAtoms = nNumAtoms;
            } else
            if ( j > nNumAtoms && (lenConnTable || j != 1) ) {
                ret = RI_ERR_SYNTAX; /* syntax error */
                goto exit_function;
            }
            /* allocate connection table */
            if ( pInChI[iComponent+i].nConnTable ) {
                inchi_free( pInChI[iComponent+i].nConnTable );
            }
            if ( !nNumAtoms && !nNumBonds && !lenConnTable ) {
                lenConnTable = 1;  /* one atom, no bonds */
            }
            pInChI[iComponent+i].nConnTable = (AT_NUMB *)inchi_calloc( lenConnTable+1, sizeof(pInChI[0].nConnTable[0]) );
            if ( !pInChI[iComponent+i].nConnTable ) {
                ret = RI_ERR_ALLOC; /* allocation failure */
                goto exit_function;
            }
            pInChI[iComponent+i].lenConnTable = lenConnTable;
        }

        /* Pass 1.2 parse a component and extract the bonds */
        num_open = 0;
        memset( last_atom, 0, sizeof(last_atom) );
        state='\0';   /* initial state */
        iBond = 0;
        p = pStart;
        pLB->len = 0;

        if ( *pbAbc == 1 ) {
            /* compressed */
            int num_neigh;
            num_open = 0;
            last_atom[num_open] = 2;
            while ( p < pEnd ) {
                if ( last_atom[num_open] > maxAtom ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                if ( isupper( UCINT *p ) ) {
                    curAtom = (AT_NUMB)inchi_strtol( p, &q, base );
                    if ( (ret = AddLinkedBond( last_atom[num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                        goto exit_function;
                    }
                    p = q;
                    if ( bPrevVersion == 1 ) {
                        while ( p < pEnd && *p == '-' ) {
                            p ++;
                            if ( (curAtom = (AT_NUMB)inchi_strtol( p, &q, base )) ) {
                                if ( (ret = AddLinkedBond( last_atom[num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                    goto exit_function;
                                }
                                p = q;
                            } else {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                        }
                    } else
                    if ( bPrevVersion == 0 && isdigit( *p ) ) {
                        num_neigh = (int)inchi_strtol( p, &q, 10 );
                        p = q;
                        while( num_neigh -- && p < pEnd ) {
                            if ( (curAtom = (AT_NUMB)inchi_strtol( p, &q, base )) ) {
                                if ( (ret = AddLinkedBond( last_atom[num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                    goto exit_function;
                                }
                                p = q;
                            } else {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                        }
                    }
                    last_atom[num_open] ++;
                } else {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
            }
        } else {
            while ( p < pEnd ) {
                /* each atom number except the first means a new bond */
                c = UCINT *p ++;
                switch ( c ) {
                case '(':
                case ')':
                case ',':
                case '-':
                    switch ( state ) {
                    case 'N':
                        state = c;
                        break;
                    default:
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    break;
                default:
                    if ( isdigit( c ) && (curAtom = (AT_NUMB)inchi_strtol( p-1, &q, 10 )) ) {
                        p = q;
                        switch( state ) {
                        case '\0':
                            last_atom[num_open] = curAtom;
                            state = 'N';
                            break;
                        case '(':
                            if ( (ret = AddLinkedBond( last_atom[num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                goto exit_function;
                            }
                            if ( ++ num_open >= LAST_AT_LEN ) {
                                ret = RI_ERR_PROGR; /* program error: buffer overflow */
                                goto exit_function;
                            }
                            last_atom[num_open] = curAtom;
                            state = 'N';
                            break;

                        case ')':
                            if ( !num_open ) {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                            if ( (ret = AddLinkedBond( last_atom[--num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                goto exit_function;
                            }
                            last_atom[num_open] = curAtom;
                            state = 'N';
                            break;

                        case ',':
                            if ( !num_open ) {
                                ret = RI_ERR_SYNTAX; /* syntax error */
                                goto exit_function;
                            }
                            if ( (ret = AddLinkedBond( last_atom[num_open-1], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                goto exit_function;
                            }
                            last_atom[num_open] = curAtom;
                            state = 'N';
                            break;
                        case '-':
                            if ( (ret = AddLinkedBond( last_atom[num_open], curAtom, (AT_NUMB)nNumAtoms, pLB )) ) {
                                goto exit_function;
                            }
                            last_atom[num_open] = curAtom;
                            state = 'N';
                            break;
                        default:
                            ret = RI_ERR_SYNTAX; /* syntax error */
                            goto exit_function;
                        }
                    } else {
                        ret = RI_ERR_SYNTAX; /* syntax error */
                        goto exit_function;
                    }
                    break;
                }
            }
        }
        /* store the bonds in connection table */
        if ( lenConnTable > 1 ) {
            for ( i = 0, m = 0; i < nNumAtoms; i ++ ) {
                k = 0;
                if ( (j = pLB->pBond[i+1].prev) ) {
                    while( k < MAXVAL ) {
                        neighbor[k++] = pLB->pBond[j].neigh;
                        if ( j == i+1 )
                            break;
                        j = pLB->pBond[j].prev;
                    }
                }
                if ( j != i+1 ) {
                    ret = RI_ERR_SYNTAX; /* syntax error */
                    goto exit_function;
                }
                /* sort the neighbors */
                insertions_sort_AT_NUMB( neighbor, k );
                pInChI[iComponent].nConnTable[m ++] = i+1; /* atom number */
                for ( j = 0; j < k && (int)neighbor[j] <= i; j ++ ) {
                    pInChI[iComponent].nConnTable[m ++] = neighbor[j];
                }
            }
            if ( m != lenConnTable ) {
                ret = RI_ERR_PROGR; /* program error */
                goto exit_function;
            }
        } else {
            pInChI[iComponent].nConnTable[0] = 1; /* single atom */
        }
        /* duplicate if needed */
        for ( i = 1; i < mpy_component; i ++ ) {
            /*
            if ( pInChI[iComponent+i].nConnTable ) {
                inchi_free( pInChI[iComponent+i].nConnTable );
            }
            pInChI[iComponent+i].nConnTable = (AT_NUMB *)inchi_calloc( lenConnTable+1, sizeof(pInChI[0].nConnTable[0]) );
            if ( !pInChI[iComponent+i].nConnTable ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            */
            if ( !pInChI[iComponent+i].nConnTable || pInChI[iComponent+i].lenConnTable != lenConnTable ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            memcpy ( pInChI[iComponent+i].nConnTable, pInChI[iComponent].nConnTable, lenConnTable*sizeof(pInChI[0].nConnTable[0]));
        }
        /* prepare for the next connection table */
        iComponent += i;
        if ( *pEnd )
            pStart = pEnd+1;
        else
            break;

    }
    ret = iComponent;

exit_function:
    if ( pLB->pBond ) {
        INCHI_HEAPCHK
        inchi_free( pLB->pBond ); 
    }
    return ret;   

#undef LAST_AT_LEN
}
/****************************************************************************************/
int nFillOutProtonMobileH( INChI *pInChI )
{
    int len = 1;
    pInChI->bDeleted = 1;
    /* formula */
    if ( !pInChI->szHillFormula && 
         !( pInChI->szHillFormula = (char *) inchi_calloc( len+1, sizeof(pInChI->szHillFormula[0]) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    strcpy( pInChI->szHillFormula, "H" );
    pInChI->nNumberOfAtoms = 1;

    /* atoms */
    if ( !pInChI->nAtom &&
         !(pInChI->nAtom = (U_CHAR *) inchi_calloc( len+1, sizeof(pInChI->nAtom[0]) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    pInChI->nAtom[0] = 1;
    /* charge */
    pInChI->nTotalCharge = 1;
    /* connection table */
    if ( !pInChI->nConnTable &&
         !(pInChI->nConnTable = (AT_NUMB *) inchi_calloc( len+1, sizeof(pInChI->nConnTable[0]) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    pInChI->nConnTable[0] = 1;
    pInChI->lenConnTable  = len;
    /* tautomer */
    if ( !pInChI->nTautomer &&
         !(pInChI->nTautomer = (AT_NUMB *) inchi_calloc( len+1, sizeof(pInChI->nTautomer[0]) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    /* nNum_H */
    if ( !pInChI->nNum_H &&
         !(pInChI->nNum_H = (S_CHAR *) inchi_calloc( len+1, sizeof(pInChI->nNum_H[0]) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    pInChI->nNum_H[0] = 0;

    pInChI->nTautomer[0] = 0;
    pInChI->lenTautomer  = 1;
    return 0;
}
/****************************************************************************************/
int nProtonCopyIsotopicInfo( INChI *pInChI_to, INChI *pInChI_from )
{
    if ( pInChI_from->nNumberOfIsotopicAtoms ) {
        if ( pInChI_to->nNumberOfIsotopicAtoms &&
             pInChI_from->nNumberOfIsotopicAtoms > pInChI_to->nNumberOfIsotopicAtoms ) {

            inchi_free( pInChI_to->IsotopicAtom );
            pInChI_to->IsotopicAtom = NULL;
            pInChI_to->nNumberOfIsotopicAtoms = 0;
        }
        if ( !pInChI_to->IsotopicAtom &&
             !(pInChI_to->IsotopicAtom = 
                (INChI_IsotopicAtom *)inchi_calloc(pInChI_from->nNumberOfIsotopicAtoms, 
                                                sizeof(pInChI_to->IsotopicAtom[0]) ) ) ) {
            return RI_ERR_ALLOC;
        }
        pInChI_to->nNumberOfIsotopicAtoms = pInChI_from->nNumberOfIsotopicAtoms;
        memcpy( pInChI_to->IsotopicAtom, pInChI_from->IsotopicAtom,
                pInChI_from->nNumberOfIsotopicAtoms * sizeof(pInChI_to->IsotopicAtom[0]) );
    } else {
        if ( pInChI_to->IsotopicAtom ) inchi_free( pInChI_to->IsotopicAtom );
        pInChI_to->IsotopicAtom = NULL;
        pInChI_to->nNumberOfIsotopicAtoms = 0;
    }
    return 0;
}
/****************************************************************************************/
int ParseSegmentFormula( const char *str, int bMobileH, INChI *pInpInChI[], int pnNumComponents[] )
{
    int i, j, mpy_component, mpy_atom, len, el_number;
    int nNumComponents = 0, iComponent, nNumAtoms, nNumAtomsAndH, iAtom, nNumH, nAltMobileH = ALT_TAUT(bMobileH);
    const char *p, *q, *e, *pStart, *pEnd;
    INChI *pInChI;
    char szEl[3];

    nNumAtoms = -999; /* impossible value */

    /* Pass 1. Count components */
    pStart = str;
    while( 1 ) {
        if ( !(pEnd = strchr( pStart, '.' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        p = pStart;
        if ( isdigit( *p ) )  {
            mpy_component = (int)inchi_strtol( p, &q, 10 );
            p = q;
        } else {
            mpy_component = 1;
        }
        if ( !mpy_component )
            break;
        if ( !isupper( UCINT *p ) ) {
            break; /* not a formula layer */
        }
        if ( pEnd == p )
            break; /* zero length formula */
        nNumComponents += mpy_component;
        if ( *pEnd )
            pStart = pEnd+1;
        else
            break;
    }
    pnNumComponents[bMobileH] = nNumComponents;

#if ( FIX_DALKE_BUGS == 1 )
    if ( nNumComponents > MAX_ATOMS ) {
        return RI_ERR_SYNTAX; /* syntax error */
    }
#endif    
    /* exit or error check */
    if ( !nNumComponents ) {
        if ( !*pStart || islower( UCINT *pStart ) ) {
            INCHI_HEAPCHK
            if ( bMobileH == TAUT_NON && 0 < ( nNumComponents = pnNumComponents[nAltMobileH]) ) {
                /* allocate InChI */
                if ( !( pInChI = (INChI *)inchi_calloc( nNumComponents, sizeof(INChI) ) ) ) {
                    return RI_ERR_ALLOC; /* alloc failure */
                }
                pInpInChI[bMobileH] = pInChI;
                pnNumComponents[bMobileH] = nNumComponents;
                for ( i = 0; i < nNumComponents; i ++ ) {
                    /* copy number of atoms */
                    len = pInpInChI[bMobileH][i].nNumberOfAtoms = pInpInChI[nAltMobileH][i].nNumberOfAtoms;
                    /* copy atoms */
                    len = (len+1)*sizeof(pInpInChI[0][0].nAtom[0]);
                    if ( pInpInChI[bMobileH][i].nAtom ) {
                        inchi_free( pInpInChI[bMobileH][i].nAtom );
                    }
                    if ( (pInpInChI[bMobileH][i].nAtom = (U_CHAR *) inchi_malloc( (len + 1) * sizeof(pInpInChI[0][0].nAtom[0]) )) ) {
                        memcpy(pInpInChI[bMobileH][i].nAtom, pInpInChI[nAltMobileH][i].nAtom, len);
                        pInpInChI[bMobileH][i].nAtom[len] = 0;
                    } else {
                        return RI_ERR_ALLOC; /* alloc failure */
                    }
                    /* copy Hill formula */
                    len = strlen( pInpInChI[nAltMobileH][i].szHillFormula)+1;
                    if ( pInpInChI[bMobileH][i].szHillFormula ) {
                        inchi_free( pInpInChI[bMobileH][i].szHillFormula );
                    }
                    if ( (pInpInChI[bMobileH][i].szHillFormula = (char *) inchi_malloc( inchi_max(len,2) )) ) {
                        memcpy( pInpInChI[bMobileH][i].szHillFormula, pInpInChI[nAltMobileH][i].szHillFormula, len);
                    } else {
                        return RI_ERR_ALLOC; /* alloc failure */
                    }
                }

            } else
            if ( bMobileH == TAUT_YES ) {
                int ret;
                /* allocate InChI */
                nNumComponents = 1;
                /* InChI */
                pnNumComponents[bMobileH] = nNumComponents;
                if ( !( pInChI = (INChI *)inchi_calloc( nNumComponents, sizeof(INChI) ) ) ) {
                    return RI_ERR_ALLOC; /* alloc failure */
                }
                pInpInChI[bMobileH] = pInChI;
                ret = nFillOutProtonMobileH( pInChI );
                if ( ret < 0 ) {
                    return ret;
                }
            }
            return 0;
        }
        return RI_ERR_SYNTAX; /* syntax error */
    }
    if ( *pEnd ) {
        return RI_ERR_SYNTAX; /* syntax error */
    }

    /* allocate InChI */
    if ( !( pInpInChI[bMobileH] = (INChI *)inchi_calloc( nNumComponents, sizeof(INChI) ) ) ) {
        return RI_ERR_ALLOC; /* alloc failure */
    }
    pInChI = pInpInChI[bMobileH];

    /* Pass 2. Count elements, save formulas and elements */
    pStart = str;
    iComponent = 0;
    while( 1 ) {
        if ( !(pEnd = strchr( pStart, '.' )) ) {
            pEnd = pStart + strlen(pStart);
        }
        p = pStart;
        if ( isdigit( UCINT *p ) )  {
            mpy_component = (int)inchi_strtol( p, &q, 10 );
            p = q;
        } else {
            mpy_component = 1;
        }
#if ( FIX_DALKE_BUGS == 1 )
        if ( iComponent + mpy_component > MAX_ATOMS ) {
            return RI_ERR_SYNTAX; /* syntax error */
        }
#endif    
        len = pEnd-p;
        for ( i = 0; i < mpy_component; i ++ ) {
            if ( pInChI[iComponent+i].szHillFormula ) {
                inchi_free( pInChI[iComponent+i].szHillFormula );
            }
            pInChI[iComponent+i].szHillFormula = (char*) inchi_malloc( inchi_max(len,1)+1 );
            memcpy( pInChI[iComponent].szHillFormula, p, len );
            pInChI[iComponent+i].szHillFormula[len] = '\0';
            if ( !i ) {
                /* Pass 2.1 Parse formula and count atoms except H */
                nNumAtoms = 0;
                nNumH     = 0;
                nNumAtomsAndH = 0;
                e = pInChI[iComponent].szHillFormula;
                while ( *e ) {
                    if ( !isupper( UCINT *e ) ) {
                        return RI_ERR_SYNTAX;
                    }
                    j = 0;
                    szEl[j ++] = *e ++;
                    if ( *e && islower( UCINT *e ) )
                        szEl[j ++] = *e ++;
                    szEl[j ++] = '\0';
                    if ( *e && isdigit( UCINT *e ) ) {
                        mpy_atom = (int)inchi_strtol( e, &q, 10 );
                        e = q;
                    } else {
                        mpy_atom = 1;
                    }
                    if ( !mpy_atom ) {
                        return RI_ERR_SYNTAX;
                    }
                    if ( szEl[0] == 'H' && !szEl[1] ) {
                        nNumH += mpy_atom;
                        continue; /* ignore H in counting number of atoms */
                    }

                    nNumAtoms += mpy_atom;
                }
#if ( FIX_DALKE_BUGS == 1 )
                if ( nNumAtoms > MAX_ATOMS ) {
                    return RI_ERR_SYNTAX; /* syntax error */
                }
#endif    
                nNumAtomsAndH = nNumAtoms? nNumAtoms : (nNumH > 0);
                pInChI[iComponent+i].nNumberOfAtoms = nNumAtomsAndH;
                if ( pInChI[iComponent+i].nAtom ) {
                    inchi_free( pInChI[iComponent+i].nAtom );
                }
                pInChI[iComponent+i].nAtom = (U_CHAR *) inchi_malloc((nNumAtomsAndH+1)*sizeof(pInChI[0].nAtom[0]));
                if ( !pInChI[iComponent+i].nAtom )
                    return RI_ERR_ALLOC; /* failed allocation */
                /* Pass 2.2 Store elements; this assumes no bridging H. Bridging H will be found in connection table, /c */
                iAtom = 0;
                if ( nNumAtoms > 0 ) { 
                    e = pInChI[iComponent+i].szHillFormula;
                    while ( *e ) {
                        if ( !isupper( UCINT *e ) ) {
                            return RI_ERR_SYNTAX;
                        }
                        j = 0;
                        szEl[j ++] = *e ++;
                        if ( *e && islower( UCINT *e ) )
                            szEl[j ++] = *e ++;
                        szEl[j ++] = '\0';
                        if ( *e && isdigit( UCINT *e ) ) {
                            mpy_atom = (int)inchi_strtol( e, &q, 10 );
                            e = q;
                        } else {
                            mpy_atom = 1;
                        }
                        if ( !mpy_atom ) {
                            return RI_ERR_SYNTAX;
                        }
                        if ( szEl[0] == 'H' && !szEl[1] )
                            continue; /* ignore H */
                        el_number = get_periodic_table_number( szEl );
                        if ( el_number == ERR_ELEM ) {
                            return RI_ERR_SYNTAX; /* wrong element */
                        }
                        while ( mpy_atom -- ) {
                            if ( iAtom >= nNumAtoms ) {
                                return RI_ERR_PROGR; /* program error */
                            }
                            pInChI[iComponent+i].nAtom[iAtom ++] = (U_CHAR)el_number;
                        }
                    }
                } else
                if ( nNumH > 0 ) {
                    pInChI[iComponent+i].nAtom[iAtom ++] = EL_NUMBER_H;
                    nNumAtoms = 1;
                }
                pInChI[iComponent+i].nAtom[iAtom] = '\0';
                if ( nNumAtoms != iAtom ) {
                    return RI_ERR_PROGR; /* program error */
                }
            } else {
                /* Copy duplicated formula */
                strcpy(pInChI[iComponent+i].szHillFormula, pInChI[iComponent].szHillFormula);
                /* Copy atoms in the duplicated formula */
                pInChI[iComponent+i].nNumberOfAtoms = nNumAtoms;
                if ( pInChI[iComponent+i].nAtom ) {
                    inchi_free( pInChI[iComponent+i].nAtom );
                }
                pInChI[iComponent+i].nAtom = (U_CHAR *) inchi_malloc(nNumAtoms+1);
                if ( !pInChI[iComponent+i].nAtom )
                    return RI_ERR_ALLOC; /* failed allocation */
                memcpy( pInChI[iComponent+i].nAtom, pInChI[iComponent].nAtom, nNumAtoms+1 );
            }
        }
        iComponent += i;
        if ( *pEnd ) {
            if ( *pEnd != '.' ) {
                return RI_ERR_SYNTAX; /* syntax error */
            }
            pStart = pEnd+1;
        } else
            break;
    }
                
    if ( iComponent != nNumComponents ) {
        return RI_ERR_PROGR; /* program error */
    }
    if ( bMobileH == TAUT_NON ) {
        /* at this point the exact number of atoms including bridging H is known from TAUT_YES */
        for ( i = 0; i < nNumComponents && i < pnNumComponents[nAltMobileH]; i ++ ) {
            if ( pInpInChI[bMobileH][i].nNumberOfAtoms < (len=pInpInChI[nAltMobileH][i].nNumberOfAtoms) ) {
                /* there are bridging H in this component */
                if ( pInpInChI[nAltMobileH][i].nAtom ) {
                    U_CHAR *nAtom = (U_CHAR *) inchi_malloc( (len+1) * sizeof(nAtom[0]) );
                    if ( !nAtom ) {
                        return RI_ERR_ALLOC;
                    }
                    memcpy( nAtom, pInpInChI[nAltMobileH][i].nAtom, len*sizeof(nAtom[0]) );
                    nAtom[ len ] = 0;
                    if ( pInpInChI[bMobileH][i].nAtom ) {
                        inchi_free( pInpInChI[bMobileH][i].nAtom );
                    }
                    pInpInChI[bMobileH][i].nAtom = nAtom;
                }
                pInpInChI[bMobileH][i].nNumberOfAtoms = len;
            }
        }
    }

    return nNumComponents+1;
}

/****************************************************************************************/
int CopySegment( INChI *pInChITo, INChI *pInChIFrom, int SegmentType, int bIsotopicTo, int bIsotopicFrom)
{
    int            ret = RI_ERR_ALLOC;
    int            len;


    if ( SegmentType==CPY_SP2 ||
         SegmentType==CPY_SP3 ||
         SegmentType==CPY_SP3_M ||
         SegmentType==CPY_SP3_S ) {

        INChI_Stereo **pstereoTo = NULL;
        INChI_Stereo *stereoFrom = bIsotopicFrom==1? pInChIFrom->StereoIsotopic :
                                   bIsotopicFrom==0? pInChIFrom->Stereo : NULL;
        if ( stereoFrom || bIsotopicFrom < 0 ) {
            if ( SegmentType==CPY_SP2 ) {
                if ( bIsotopicFrom < 0 ||
                     (stereoFrom->b_parity &&
                     stereoFrom->nBondAtom1 &&
                     stereoFrom->nBondAtom2) ) {

                    len = (bIsotopicFrom < 0)? 0 : stereoFrom->nNumberOfStereoBonds;
                    pstereoTo = bIsotopicTo? &pInChITo->StereoIsotopic : &pInChITo->Stereo;
                    if ( !pstereoTo[0] ) {
                        if ( !(pstereoTo[0] = (INChI_Stereo *)inchi_calloc( 1, sizeof(**pstereoTo))) ) {
                            goto exit_function;
                        }
                    }
                    if ( pstereoTo[0]->nNumberOfStereoBonds > 0 || pstereoTo[0]->b_parity ||
                         pstereoTo[0]->nBondAtom1 || pstereoTo[0]->nBondAtom2 ) {
                        ret = RI_ERR_SYNTAX; /* stereo already exists */
                        goto exit_function;
                    }
                    /* allocate sp2 stereo */
                    if ( !(pstereoTo[0]->b_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pstereoTo[0]->b_parity[0]) ) ) ||
                         !(pstereoTo[0]->nBondAtom1 = (AT_NUMB *)inchi_calloc( len+1, sizeof(pstereoTo[0]->nBondAtom1[0]) ) ) ||
                         !(pstereoTo[0]->nBondAtom2 = (AT_NUMB *)inchi_calloc( len+1, sizeof(pstereoTo[0]->nBondAtom2[0]) ) ) ) {
                        /* cleanup */
                        if ( pstereoTo[0]->b_parity ) {
                            INCHI_HEAPCHK
                            inchi_free( pstereoTo[0]->b_parity );
                            pstereoTo[0]->b_parity = NULL;
                        }
                        if ( pstereoTo[0]->nBondAtom1 ) {
                            INCHI_HEAPCHK
                            inchi_free( pstereoTo[0]->nBondAtom1 );
                            pstereoTo[0]->nBondAtom1 = NULL;
                        }
                        if ( pstereoTo[0]->nBondAtom2 ) {
                            INCHI_HEAPCHK
                            inchi_free( pstereoTo[0]->nBondAtom2 );
                            pstereoTo[0]->nBondAtom2 = NULL;
                        }
                        INCHI_HEAPCHK
                        goto exit_function;
                    }
                    /* copy stereo */
                    if ( bIsotopicFrom >= 0 && len ) {
                        memcpy( pstereoTo[0]->b_parity, stereoFrom->b_parity, (len+1)*sizeof(pstereoTo[0]->b_parity[0]) );
                        memcpy( pstereoTo[0]->nBondAtom1, stereoFrom->nBondAtom1, (len+1)*sizeof(pstereoTo[0]->nBondAtom1[0]) );
                        memcpy( pstereoTo[0]->nBondAtom2, stereoFrom->nBondAtom2, (len+1)*sizeof(pstereoTo[0]->nBondAtom2[0]) );
                    }
                    pstereoTo[0]->nNumberOfStereoBonds = len;
                
                    return len+1;
                } else {
                    return 0;
                }
            } else
            if ( SegmentType==CPY_SP3 ) {
                if ( bIsotopicFrom < 0 ||
                     (stereoFrom->t_parity &&
                     stereoFrom->nNumber) ) {

                    len = (bIsotopicFrom < 0)? 0 : stereoFrom->nNumberOfStereoCenters;
                    pstereoTo = bIsotopicTo? &pInChITo->StereoIsotopic : &pInChITo->Stereo;
                    if ( !pstereoTo[0] ) {
                        if ( !(pstereoTo[0] = (INChI_Stereo *)inchi_calloc( 1, sizeof(**pstereoTo))) ) {
                            goto exit_function;
                        }
                    }
                    if ( pstereoTo[0]->nNumberOfStereoCenters > 0 || pstereoTo[0]->t_parity ||
                         pstereoTo[0]->nNumber ) {
                        ret = RI_ERR_SYNTAX; /* stereo already exists */
                        goto exit_function;
                    }
                    /* allocate sp3 stereo */
                    if ( !(pstereoTo[0]->t_parity = (S_CHAR *)inchi_calloc( len+1, sizeof(pstereoTo[0]->b_parity[0]) ) ) ||
                         !(pstereoTo[0]->nNumber = (AT_NUMB *)inchi_calloc( len+1, sizeof(pstereoTo[0]->nBondAtom1[0]) ) ) ) {
                        /* cleanup */
                        if ( pstereoTo[0]->t_parity ) {
                            inchi_free( pstereoTo[0]->t_parity );
                            pstereoTo[0]->t_parity = NULL;
                        }
                        if ( pstereoTo[0]->nNumber ) {
                            inchi_free( pstereoTo[0]->nNumber );
                            pstereoTo[0]->nNumber = NULL;
                        }
                        goto exit_function;
                    }
                    /* copy stereo */
                    if ( bIsotopicFrom >= 0 && len ) {
                    memcpy( pstereoTo[0]->t_parity, stereoFrom->t_parity, (len+1)*sizeof(pstereoTo[0]->t_parity[0]) );
                    memcpy( pstereoTo[0]->nNumber, stereoFrom->nNumber, (len+1)*sizeof(pstereoTo[0]->nNumber[0]) );
                    }
                    pstereoTo[0]->nNumberOfStereoCenters = len;
                    return len+1;
                } else {
                    return 0;
                }
            } else
            if ( SegmentType==CPY_SP3_M ) {
                pstereoTo = bIsotopicTo? &pInChITo->StereoIsotopic : &pInChITo->Stereo;
                if ( !pstereoTo[0] ) {
                    if ( !(pstereoTo[0] = (INChI_Stereo *)inchi_calloc( 1, sizeof(**pstereoTo))) ) {
                        goto exit_function;
                    }
                }
                if ( pstereoTo[0]->nCompInv2Abs && NO_VALUE_INT != pstereoTo[0]->nCompInv2Abs ) {
                    ret = RI_ERR_SYNTAX; /* stereo already exists */
                    goto exit_function;
                }
                if ( bIsotopicFrom < 0 ) {
                    pstereoTo[0]->nCompInv2Abs = 0;
                } else {
                    pstereoTo[0]->nCompInv2Abs = stereoFrom->nCompInv2Abs;
                }
                return 1;
            } else
            /* use bTrivialInv to save /s1, /s2, /s3 */
            if ( SegmentType==CPY_SP3_S ) {
                pstereoTo = bIsotopicFrom? &pInChITo->StereoIsotopic : &pInChITo->Stereo;
                if ( !pstereoTo[0] ) {
                    if ( !(pstereoTo[0] = (INChI_Stereo *)inchi_calloc( 1, sizeof(**pstereoTo))) ) {
                        goto exit_function;
                    }
                }
                if ( pstereoTo[0]->bTrivialInv ) {
                    ret = RI_ERR_SYNTAX; /* stereo already exists */
                    goto exit_function;
                }
                pstereoTo[0]->bTrivialInv = stereoFrom->bTrivialInv;
                if ( bIsotopicFrom < 0 ) {
                    pstereoTo[0]->bTrivialInv = 0;
                } else {
                    pstereoTo[0]->bTrivialInv =  stereoFrom->bTrivialInv;
                }
                return 1;
            }
        }
        return 0; /* nothing to copy */
    } else
    if ( SegmentType == CPY_ISO_AT ) {
        int                  nNumberOfIsotopicAtoms = pInChIFrom->nNumberOfIsotopicAtoms;
        INChI_IsotopicAtom   **pIsotopicAtomTo = NULL;
        INChI_IsotopicAtom    *IsotopicAtomFrom = pInChIFrom->IsotopicAtom;
        if ( bIsotopicFrom < 0 || IsotopicAtomFrom ) {
            len = (bIsotopicFrom < 0)? 0 : nNumberOfIsotopicAtoms;
            pIsotopicAtomTo = &pInChITo->IsotopicAtom;
            if ( !*pIsotopicAtomTo ) {
                if ( !( *pIsotopicAtomTo = (INChI_IsotopicAtom *)inchi_calloc( len+1, sizeof(**pIsotopicAtomTo) ) ) ) {
                    goto exit_function;
                }
            }
            if ( pInChITo->nNumberOfIsotopicAtoms ) {
                ret = RI_ERR_SYNTAX; /* stereo already exists */
                goto exit_function;
            }
            if ( bIsotopicFrom >= 0 && len ) {
                memcpy( *pIsotopicAtomTo, IsotopicAtomFrom, (len+1)*sizeof(**pIsotopicAtomTo) );
            }
            pInChITo->nNumberOfIsotopicAtoms = len;
            return len+1;
        }
        return 0;
    }
    ret = RI_ERR_PROGR; /* program error */
exit_function:
    return ret;


}

/**********************************************************************************/
/*  Sort neighbors in ascending order */
int insertions_sort_AT_NUMB( AT_NUMB *base, int num )
{
  AT_NUMB *i, *j, *pk, tmp;
  int  k, num_trans = 0;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1, tmp = *j; j > base && *i > tmp; j=i, i -- ) {
        *j = *i;
        num_trans ++;
     }
     *j = tmp;
  }
  return num_trans;
}


/* read */

int getInChIChar( INCHI_IOSTREAM *pInp )
{
    if (pInp->type==INCHI_IOSTREAM_STRING)
    {
        /* input from string */
        if ( pInp->s.nPtr < pInp->s.nUsedLength ) 
            return (int) pInp->s.pStr[pInp->s.nPtr++];
        return RI_ERR_EOF;
    }

    else
    {
        /* input from plain file */
        int c;
#if ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        do 
        {
            c = getc( pInp->f );
            if ( c == EOF ) 
            {
                c = RI_ERR_EOF;
                break;
            }
        } 
        while( c == '\r' );
#else
        c = getc( pInp->f );
        if ( c == EOF ) 
        {
            c = RI_ERR_EOF;
        }
#endif
        return c;
    }

    
}

int AddInChIChar( INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken )
{
    int c = getInChIChar( pInp );
    /*
    while ( c == '\r' ) {
        c = getInChIChar( pInp );
    }
    */
    INCHI_HEAPCHK
    if ( Line->len + 2 >= Line->len_alloc ) {
        char *str = (char *) inchi_calloc( Line->len_alloc + SEGM_LINE_ADD, sizeof(str[0]) );
        INCHI_HEAPCHK
        if ( str ) {
            if ( Line->len > 0 && Line->str ) {
                memcpy( str, Line->str, sizeof(str[0]) * Line->len );
                Line->len_alloc += SEGM_LINE_ADD;
                inchi_free( Line->str );
                INCHI_HEAPCHK
            } else {
                Line->len_alloc += SEGM_LINE_ADD;
            }
            Line->str = str;
        } else {
            c = RI_ERR_ALLOC; /* fatal error */
            goto exit_function;
        }
    }
    INCHI_HEAPCHK
    if ( c < 0 ) {
        Line->str[Line->len] = '\0';
        INCHI_HEAPCHK
        c = RI_ERR_SYNTAX; /* fatal error: wrong char */
        goto exit_function;
    }
    if ( c && strchr( pszToken, c ) ) {
        Line->str[Line->len] = '\0';
        INCHI_HEAPCHK
        c = -(c+2);
        goto exit_function;
    } else
    if ( !c && !Line->len ) {
        Line->str[Line->len] = c;
        INCHI_HEAPCHK
    } else {
        Line->str[Line->len ++] = c;
        INCHI_HEAPCHK
    }
exit_function:
    INCHI_HEAPCHK
    return c;
}
int nGetInChISegment( INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken )
{
    int c;
    Line->len = 0;
    while( 0 < (c = AddInChIChar( pInp, Line, pszToken ) ) )
        ;
    if ( c < - 2 ) {
        c = -(c+2);
    }
    Line->c = c;
    return c;
}
/********************************************************************************/
/* add one more bond to the linked lists for both neighbors */
int AddLinkedBond( AT_NUMB at1, AT_NUMB at2, AT_NUMB num_at, LINKED_BONDS *pLB )
{
    int nReqLen = inchi_max( 2*num_at+2, pLB->len + 2 );
    AT_NUMB prev;
    if ( pLB->len_alloc <= nReqLen ) {
        /*int nNewLen = nReqLen + (nReqLen + LINKED_BOND_ADD - 1)%LINKED_BOND_ADD + LINKED_BOND_ADD;*/
        int nNewLen = nReqLen - nReqLen%LINKED_BOND_ADD + 2*LINKED_BOND_ADD;
        ONE_LINKED_BOND *pBond = (ONE_LINKED_BOND *)inchi_calloc( nNewLen, sizeof(pBond[0]) );
        if ( !pBond )
            return RI_ERR_ALLOC; /* allocation error */
        if ( pLB->pBond && pLB->len ) {
            memcpy( pBond, pLB->pBond, pLB->len*sizeof(pBond[0]) );
        }
        if ( pLB->pBond )
            inchi_free( pLB->pBond );
        pLB->pBond = pBond;
        pLB->len_alloc = nNewLen;
    }
    if ( !pLB->len ) {
        pLB->len = num_at+1;
        memset( pLB->pBond, 0, (num_at+1)*sizeof(pLB->pBond[0]) );
    }
    
    prev = pLB->pBond[at1].prev; /* position of the last neighbor of at1 in the pLB->pBond */
    if ( !prev ) {
        pLB->pBond[at1].neigh = at2;
        pLB->pBond[at1].prev  = at1;
    } else {
        pLB->pBond[pLB->len].neigh = at2;
        pLB->pBond[pLB->len].prev  = prev;
        pLB->pBond[at1].prev       = pLB->len ++;
    }
    
    prev = pLB->pBond[at2].prev; /* position of the last neighbor of at2 in the pLB->pBond */
    if ( !prev ) {
        pLB->pBond[at2].neigh = at1;
        pLB->pBond[at2].prev  = at2;
    } else {
        pLB->pBond[pLB->len].neigh = at1;
        pLB->pBond[pLB->len].prev  = prev;
        pLB->pBond[at2].prev       = pLB->len ++;
    }
    return 0;
}

#endif /* READ_INCHI_STRING */
