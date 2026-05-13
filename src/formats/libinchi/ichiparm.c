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

 /* djb-rwth: used for strdup function */
#ifndef _WIN32
#if defined(__STDC_ALLOC_LIB__)
#define __STDC_WANT_LIB_EXT2__ 1
#else
#define _POSIX_C_SOURCE 200809L
#endif
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
 /* #include <varargs.h> */
#include <errno.h>
#include <limits.h>

#include "mode.h"
#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif

#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"
#ifdef TARGET_LIB_FOR_WINCHI
#include "../../../IChI_lib/src/ichi_lib.h"
#endif
#include "ichicomp.h"

#include "ichi_io.h"
#include "util.h"

#include "bcf_s.h"

#define VER103_DEFAULT_MODE    (REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO |\
                                REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)


int DetectInputINChIFileType(FILE** inp_file, INPUT_PARMS* ip, const char* fmode);

int set_common_options_by_parg(const char* pArg,
    int  developer_options,
    INPUT_PARMS* ip,
    INCHI_MODE* pbVer1DefaultMode,
    int* pnMode,
    int* pbINChIOutputOptions,
    int* pbINChIOutputOptions2,
    int* pbStdFormat,
    int* pbHashKey,
    int* pbHashXtra1,
    int* pbHashXtra2,
    int* pbFixSp3bug,
    int* pbFixFB2,
    int* pbAddPhosphineStereo,
    int* pbAddArsineStereo,
    int* pbNoStructLabels,
    int* pbPointedEdgeStereo,
    int* pbDoNotAddH,
    int* pbForcedChiralFlag,
    int* pbReconnectCoord,
    int* pbKetoEnolTaut,
    int* pb15TautNonRing,
    int* pbPT_06_00_Taut,
    int* pbPT_13_00_Taut,
    int* pbPT_16_00_Taut,
    int* pbPT_18_00_Taut,
    int* pbPT_22_00_Taut,
    int* pbPT_39_00_Taut,
    int* pbLooseTSACheck,
    int* pbLargeMolecules,
    int* pbPolymers,
    int* pbFoldPolymerSRU,
    int* pbFrameShiftScheme,
    int* pbStereoAtZz,
    int* pbNPZz,
    int* pbNoWarnings,
    int* pbMergeHash,
    int* pbHideInChI);


/****************************************************************************
Returns 1 if pArg recognized and treated, 0 otherwise
****************************************************************************/
int set_common_options_by_parg(const char* pArg,
    int  developer_options,
    INPUT_PARMS* ip,
    INCHI_MODE* pbVer1DefaultMode,
    int* pnMode,
    int* pbINChIOutputOptions,
    int* pbINChIOutputOptions2,
    int* pbStdFormat,
    int* pbHashKey,
    int* pbHashXtra1,
    int* pbHashXtra2,
    int* pbFixSp3bug,
    int* pbFixFB2,
    int* pbAddPhosphineStereo,
    int* pbAddArsineStereo,
    int* pbNoStructLabels,
    int* pbPointedEdgeStereo,
    int* pbDoNotAddH,
    int* pbForcedChiralFlag,
    int* pbReconnectCoord,
    int* pbKetoEnolTaut,
    int* pb15TautNonRing,
    int* pbPT_06_00_Taut,
    int* pbPT_13_00_Taut,
    int* pbPT_16_00_Taut,
    int* pbPT_18_00_Taut,
    int* pbPT_22_00_Taut,
    int* pbPT_39_00_Taut,
    int* pbLooseTSACheck,
    int* pbLargeMolecules,
    int* pbPolymers,
    int* pbFoldPolymerSRU,
    int* pbFrameShiftScheme,
    int* pbStereoAtZz,
    int* pbNPZz,
    int* pbNoWarnings,
    int* pbMergeHash,
    int* pbHideInChI
)
{
    int got = 0;

    /* Input */
    if (!inchi_stricmp(pArg, "INPAUX"))
    {
        if (INPUT_NONE == ip->nInputType)
        {
            ip->nInputType = INPUT_INCHI_PLAIN;
        }
        got = 1;
    }
    else if (!inchi_memicmp(pArg, "START:", 6))
    {
        ip->first_struct_number = strtol(pArg + 6, NULL, 10);
        got = 1;
    }
    else if (!inchi_memicmp(pArg, "END:", 4))
    {
        ip->last_struct_number = strtol(pArg + 4, NULL, 10);
        got = 1;
    }
    else if (!inchi_memicmp(pArg, "RECORD:", 7))
    {
        long num = strtol(pArg + 7, NULL, 10);
        /* djb-rwth: removing redundant code */
        ip->first_struct_number = num;
        ip->last_struct_number = num;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NOLABELS"))
    {
        *pbNoStructLabels = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "SAVEOPT"))
    {
        (*pbINChIOutputOptions) |= INCHI_OUT_SAVEOPT;
        got = 1;
    }
    /* Generation */
    else if (!inchi_stricmp(pArg, "AUXNONE"))
    {
        /* no aux. info */
        (*pbINChIOutputOptions) |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
        (*pbINChIOutputOptions) &= ~INCHI_OUT_SHORT_AUX_INFO;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "MISMATCHISERROR"))
    {
        /* Consider InChI conversion "problem/mismatch" as error */
        (*pbINChIOutputOptions2) |= INCHI_OUT_MISMATCH_AS_ERROR;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "OUTERRINCHI"))
    {
        /* Signify InChI error generation on InChI strings output, not only report to log file */
        (*pbINChIOutputOptions2) |= INCHI_OUT_INCHI_GEN_ERROR;
        got = 1;
    }
    /* InChIKey/InChI hash */
    else if (!inchi_stricmp(pArg, "Key"))
    {
        *pbHashKey = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "XHash1"))
    {
        *pbHashXtra1 = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "XHash2"))
    {
        *pbHashXtra2 = 1;
        got = 1;
    }
    /* All modes (std and non-std InChI) structure perception options */
    /* These options DO NOT TURN OFF Std flag                         */
    else if (!inchi_stricmp(pArg, "SNON"))
    {
        (*pbVer1DefaultMode) &= ~REQ_MODE_STEREO; /* no stereo */
        (*pnMode) &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NEWPSOFF"))
    {
        *pbPointedEdgeStereo = 0;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "DONOTADDH"))
    {
        *pbDoNotAddH = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "LooseTSACheck"))
    {
        (*pbLooseTSACheck) = 1;
        got = 1;
    }

#ifndef USE_STDINCHI_API
    /* These options DO TURN OFF Std flag   */
    else if (!inchi_stricmp(pArg, "SREL"))
    {
        if ((*pnMode) & REQ_MODE_RACEMIC_STEREO)
        {
            (*pnMode) ^= REQ_MODE_RACEMIC_STEREO;
        }
        if ((*pnMode) & REQ_MODE_CHIR_FLG_STEREO)
        {
            (*pnMode) ^= REQ_MODE_CHIR_FLG_STEREO;
        }
        (*pnMode) |= REQ_MODE_RELATIVE_STEREO;
        (*pnMode) |= REQ_MODE_STEREO;
        *pbStdFormat = 0;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "SRAC"))
    {
        if ((*pnMode) & REQ_MODE_RELATIVE_STEREO)
        {
            (*pnMode) ^= REQ_MODE_RELATIVE_STEREO;
        }
        if ((*pnMode) & REQ_MODE_CHIR_FLG_STEREO)
        {
            (*pnMode) ^= REQ_MODE_CHIR_FLG_STEREO;
        }
        (*pnMode) |= REQ_MODE_RACEMIC_STEREO;
        (*pnMode) |= REQ_MODE_STEREO;
        *pbStdFormat = 0;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "SUCF"))
    {
        if ((*pnMode) & REQ_MODE_RELATIVE_STEREO)
        {
            (*pnMode) ^= REQ_MODE_RELATIVE_STEREO;
        }
        if ((*pnMode) & REQ_MODE_RACEMIC_STEREO)
        {
            (*pnMode) ^= REQ_MODE_RACEMIC_STEREO;
        }
        (*pnMode) |= REQ_MODE_CHIR_FLG_STEREO; /* stereo defined by the Chiral flag */
        (*pnMode) |= REQ_MODE_STEREO;
        *pbStdFormat = 0;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "ChiralFlagON"))
    {
        /* used only with /SUCF */
        /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
        (*pbForcedChiralFlag) &= ~FLAG_SET_INP_AT_NONCHIRAL;
        (*pbForcedChiralFlag) |= FLAG_SET_INP_AT_CHIRAL;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "ChiralFlagOFF"))
    {
        /* used only with /SUCF */
        /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
        (*pbForcedChiralFlag) &= ~FLAG_SET_INP_AT_CHIRAL;
        (*pbForcedChiralFlag) |= FLAG_SET_INP_AT_NONCHIRAL;
        got = 1;
    }

    /*--- Non-std InChI creation options ---*/
    /* These options DO TURN OFF Std flag       */
    else if (!inchi_stricmp(pArg, "SUU"))
    {
        /* include omitted undef/unknown stereo */
        (*pbVer1DefaultMode) &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
        *pbStdFormat = 0;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "SLUUD"))
    {
        /* make labels for unknown and undefined stereo different */
        (*pbVer1DefaultMode) |= REQ_MODE_DIFF_UU_STEREO;
        *pbStdFormat = 0;
        got = 1;
    }
    /* FixedH */
    else if (!inchi_stricmp(pArg, "FIXEDH"))
    {
        (*pbVer1DefaultMode) |= REQ_MODE_BASIC;  /* tautomeric */
        *pbStdFormat = 0;
        got = 1;
    }
    /* RecMet */
    else if (!inchi_stricmp(pArg, "RECMET"))
    {
        /* reconnect metals */
        *pbReconnectCoord = 1;
        *pbStdFormat = 0;
        got = 1;
    }
#if ( KETO_ENOL_TAUT == 1 )
    else if (!inchi_stricmp(pArg, "KET"))
    {
        *pbKetoEnolTaut = 1;
        *pbStdFormat = 0;
        got = 1;
    }
#endif
#if ( TAUT_15_NON_RING == 1 )
    else if (!inchi_stricmp(pArg, "15T"))
    {
        *pb15TautNonRing = 1;
        *pbStdFormat = 0;
        got = 1;
    }
#endif

#if ( TAUT_PT_22_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_22_00"))
    {
        *pbPT_22_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

#if ( TAUT_PT_16_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_16_00"))
    {
        *pbPT_16_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

#if ( TAUT_PT_06_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_06_00"))
    {
        *pbPT_06_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

#if ( TAUT_PT_39_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_39_00"))
    {
        *pbPT_39_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

#if ( TAUT_PT_13_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_13_00"))
    {
        *pbPT_13_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

#if ( TAUT_PT_18_00  == 1 )
    else if (!inchi_stricmp(pArg, "PT_18_00"))
    {
        *pbPT_18_00_Taut = 1;
        *pbStdFormat = 0;
    }
#endif

    else if (!inchi_stricmp(pArg, "LargeMolecules"))
    {
        *pbLargeMolecules = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "Polymers"))
    {
        *pbPolymers = POLYMERS_MODERN;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "Polymers105"))
    {
        *pbPolymers = POLYMERS_LEGACY;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NPZz"))
    {
        *pbNPZz = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NoWarnings"))
    {
        *pbNoWarnings = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "MergeHash"))
    {
        *pbMergeHash = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NoInChI") || !inchi_stricmp(pArg, "HideInChI"))
    {
        *pbHideInChI = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "FoldCRU")) /* v. 1.06 */
    {
        *pbFoldPolymerSRU = 1;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "FoldSRU")) /* v. 1.06 */
    {
        *pbFoldPolymerSRU = 1;
        got = 1;
    }
    else if (!inchi_memicmp(pArg, "FrameShift:", 11))
    {
        int k;
        char wrd[256];
        k = 0;
        mystrncpy(wrd, pArg + 11, 256);
        lrtrim(wrd, &k);
        if (k)
        {
            if (!inchi_stricmp(wrd, "None"))
            {
                *pbFrameShiftScheme = FSS_NONE;
            }
            else if (!inchi_stricmp(wrd, "Cyclize"))
            {
                *pbFrameShiftScheme = FSS_STARS_CYCLED;
            }
            else if (!inchi_stricmp(wrd, "MoveStars"))
            {
                *pbFrameShiftScheme = FSS_STARS_OPENED;
            }
            else if (!inchi_stricmp(wrd, "MoveBrackets"))
            {
                *pbFrameShiftScheme = FSS_STARS_ENDS_OPENED;
            }
        }
        else
        {
            *pbFrameShiftScheme = FSS_STARS_CYCLED;
        }
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NoFrameShift"))
    {
        *pbFrameShiftScheme = FSS_NONE;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "NoEdits"))
    {
        *pbFoldPolymerSRU = 0;
        *pbFrameShiftScheme = FSS_NONE;
        got = 1;
    }
    else if (!inchi_stricmp(pArg, "SATZZ"))
    {
        *pbStereoAtZz = 1;
        got = 1;
    }


#endif /* ifndef USE_STDINCHI_API */

    if (!got && developer_options)
    {

        if (!inchi_stricmp(pArg, "PGO"))
        {
            /* PGO : extract all good MOLfiles into the problem file */
            ip->bSaveAllGoodStructsAsProblem = 1;
            got = 1;
        }
#if ( ALLOW_SUBSTRUCTURE_FILTERING== 1 )
        else if (!inchi_stricmp(pArg, "FilterSS"))
        {
            ip->bFilterSS = 1;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "InvFilterSS"))
        {
            ip->bFilterSS = -1;
            got = 1;
        }
#endif
        /* Options below DO TURN OFF Std flag   */
        if (!inchi_stricmp(pArg, "FNUDOFF"))
        {
            ip->bFixNonUniformDraw = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "FixSp3bugOFF"))
        {
            *pbFixSp3bug = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "FBOFF"))
        {
            *pbFixSp3bug = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "FB2OFF"))
        {
            *pbFixFB2 = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "SPXYZOFF"))
        {
            *pbAddPhosphineStereo = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "SASXYZOFF"))
        {
            *pbAddArsineStereo = 0;
            *pbStdFormat = 0;
            got = 1;
        }
        else if (!inchi_stricmp(pArg, "Polymers105+"))
        {
            *pbPolymers = POLYMERS_LEGACY_PLUS;
            *pbStdFormat = 0;
            got = 1;
        }


    }

    return got;
}


/****************************************************************************
  Parse command line and set processing options
****************************************************************************/
int ReadCommandLineParms(int argc,
    const char* argv[],
    INPUT_PARMS* ip,
    char* szSdfDataValue,
    unsigned long* ulDisplTime,
    int bReleaseVersion,
    INCHI_IOSTREAM* log_file)
{
#if (BUILD_WITH_ENG_OPTIONS==1)
    const int developer_options = 1;
#else
    const int developer_options = 0;
#endif
    const char* q;
    const char* ext[MAX_NUM_PATHS + 1];
    const char* pArg;
    char szNameSuffix[32], szOutputPath[512];

    unsigned long ul;

    double t = 0;

    int bVer1Options = 1;
    int nMode = 0;
    int nReleaseMode = nMode | (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO);
    INCHI_MODE bVer1DefaultMode = VER103_DEFAULT_MODE;

    int bNameSuffix;
    int bOutputPath;
    int bMergeAllInputStructures;
    int bForcedChiralFlag = 0;

    int bDisconnectSalts = (DISCONNECT_SALTS == 1);
    int bDoNotAddH = 0;
    int bRecognizedOption;
    int bTgFlagVariableProtons = 1;
    int bTgFlagHardAddRenProtons = 1;
    int bReconnectCoord = (RECONNECT_METALS == 1);
    int bDisconnectCoord = (DISCONNECT_METALS == 1);
    int bDisconnectCoordChkVal = (CHECK_METAL_VALENCE == 1);
    int bMovePositiveCharges = (MOVE_CHARGES == 1);
    int bAcidTautomerism = (DISCONNECT_SALTS == 1) ? (TEST_REMOVE_S_ATOMS == 1 ? 2 : 1) : 0;
    int bUnchargedAcidTaut = (CHARGED_SALTS_ONLY == 0);
    int bMergeSaltTGroups = (DISCONNECT_SALTS == 1);
#if ( MIN_SB_RING_SIZE > 0 )
    int nMinDbRinSize = MIN_SB_RING_SIZE, mdbr = 0;
#endif
#ifdef STEREO_WEDGE_ONLY
    int bPointedEdgeStereo = STEREO_WEDGE_ONLY; /* NEWPS TG_FLAG_POINTED_EDGE_STEREO */
#endif
#if ( FIX_ADJ_RAD == 1 )
    int bFixAdjacentRad = 0;
#endif
    int bAddPhosphineStereo = 1;
    int bAddArsineStereo = 1;
    int bFixSp3bug = 1;
    int bFixFB2 = 1;
    int bKetoEnolTaut = 0;
    int b15TautNonRing = 0;
    int bPT_22_00_Taut = 0;
    int bPT_16_00_Taut = 0;
    int bPT_06_00_Taut = 0;
    int bPT_39_00_Taut = 0;
    int bPT_13_00_Taut = 0;
    int bPT_18_00_Taut = 0;
    int bStdFormat = 1;
    int bHashKey = 0;
    int bHashXtra1 = 0;
    int bHashXtra2 = 0;
    int bLargeMolecules = 0;
    int bPolymers = POLYMERS_NO;
#ifdef TARGET_LIB_FOR_WINCHI
    int bFoldPolymerSRU = 0;
    int bFrameShiftScheme = FSS_STARS_CYCLED;
#else
    int bFoldPolymerSRU = 0;
    int bFrameShiftScheme = FSS_STARS_CYCLED;
#endif
    int bLooseTSACheck = 0;
    int bStereoAtZz = 0;
    int bNPZz = 0;
    int bNoWarnings = 0;
    int bMergeHash = 0;
    int bHideInChI = 0;

    int bOutputStyle = INCHI_OUT_PLAIN_TEXT;
    int bDisplay = 0;
    int bNoStructLabels = 0;
    int bOutputMolfileOnly = 0;
    int bOutputMolfileDT = 0;
    int bOutputMolfileSplit = 0;
    int bDisplayCompositeResults = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    int nFontSize = -9; /* djb-rwth: ignoring LLVM warning: variable used */
    int bINChIOutputOptions2 = 0;
#ifdef TARGET_LIB_FOR_WINCHI
    int is_gui = 1;
    int bINChIOutputOptions = INCHI_OUT_EMBED_REC; /* embed reconnected & output full aux info */
    int bCompareComponents = CMP_COMPONENTS;
#else
    int is_gui = 0;
    int bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0);
    int bCompareComponents = 0;
#endif
#if ( READ_INCHI_STRING == 1 )
    int bDisplayIfRestoreWarnings = 0;
#endif
#if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
    int numOutNameExt;
    char szOutNameExt[3][128];
#endif

    int i, k, c, got;
    int timeout_set_warning = 0;
    int timeout_set_error = 0;

    ext[0] = ".mol";
    ext[1] = bVer1Options ? ".txt" : ".ich";
    ext[2] = ".log";
    ext[3] = ".prb";
    ext[4] = "";

    /*  Init table of parameters */
    memset(ip, 0, sizeof(*ip)); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Default are StdInChI generation options */
    bVer1DefaultMode &= ~REQ_MODE_BASIC; /* "FIXEDH - OFF" */
    bReconnectCoord = 0;                /* "RECMET - OFF" */
    bPointedEdgeStereo = 1;                /* "NEWPS"        */
    ip->bFixNonUniformDraw = 1;                /* "FNUD"         */

#ifndef COMPILE_ANSI_ONLY
    strcpy(ip->tdp.ReqShownFoundTxt[ilSHOWN], "Shown");
    ip->dp.sdp.tdp = &ip->tdp;
    ip->dp.pdp = &ip->pdp;
#endif
    memset(szNameSuffix, 0, sizeof(szNameSuffix)); /* djb-rwth: memset_s C11/Annex K variant? */
    bNameSuffix = 0;
    memset(szOutputPath, 0, sizeof(szOutputPath)); /* djb-rwth: memset_s C11/Annex K variant? */
    bOutputPath = 0;
#if( OUTPUT_FILE_EXT == 1 )
    memset(szOutNameExt, 0, sizeof(szOutNameExt)); /* djb-rwth: memset_s C11/Annex K variant? */
    numOutNameExt = 0;
#endif
    bMergeAllInputStructures = 0;
#ifdef TARGET_API_LIB
    ip->msec_MaxTime = 0;      /*  milliseconds, default in libinchi: unlimited */
#else
    ip->msec_MaxTime = 60000;  /*  milliseconds, default: 60 sec */
#endif
    * ulDisplTime = 0;

    if (bReleaseVersion)
    {
        /*  normal */
        ip->bAbcNumbers = 0;
        ip->bCtPredecessors = 0;
    }
    else
    {
        nReleaseMode = 0;
    }
    if (bVer1Options)
    {
        bNameSuffix = 1;
        szNameSuffix[0] = '\0';
    }

#if ( ALLOW_SUBSTRUCTURE_FILTERING== 1 )
    ip->bFilterSS = 0;
#endif



    /* Analyze command line switches */
    for (i = 1; i < argc; i++)
    {

        if (is_gui && INCHI_OPTION_PREFX == argv[i][0] && INCHI_OPTION_PREFX != argv[i][1])
        {
            /* Parsing TARGET_LIB_FOR_WINCHI GUI options (and v. 0.9xx Beta as well)  */
            pArg = argv[i] + 1;
            got = set_common_options_by_parg(pArg, developer_options, ip, &bVer1DefaultMode, &nMode,
                &bINChIOutputOptions, &bINChIOutputOptions2,
                &bStdFormat, &bHashKey, &bHashXtra1, &bHashXtra2,
                &bFixSp3bug, &bFixFB2,
                &bAddPhosphineStereo, &bAddArsineStereo,
                &bNoStructLabels, &bPointedEdgeStereo,
                &bDoNotAddH, &bForcedChiralFlag, &bReconnectCoord,
                &bKetoEnolTaut,
                &b15TautNonRing,
                &bPT_06_00_Taut, &bPT_13_00_Taut, &bPT_16_00_Taut,
                &bPT_18_00_Taut, &bPT_22_00_Taut, &bPT_39_00_Taut,
                &bLooseTSACheck,
                &bLargeMolecules, &bPolymers,
                &bFoldPolymerSRU, &bFrameShiftScheme,
                &bStereoAtZz, &bNPZz,
                &bNoWarnings, &bMergeHash, &bHideInChI);
            if (got)
            {
                ;
            }
            else if (INPUT_NONE == ip->nInputType &&
                (!inchi_memicmp(pArg, "SDF", 3)) &&
                (pArg[3] == ':'))
            {
                k = 0;
                mystrncpy(ip->szSdfDataHeader, pArg + 4, MAX_SDF_HEADER + 1);
                lrtrim(ip->szSdfDataHeader, &k);
                if (k)
                {
                    ip->pSdfLabel = ip->szSdfDataHeader;
                    ip->pSdfValue = szSdfDataValue;
                    ip->nInputType = INPUT_SDFILE;
                }
                else
                {
                    ip->pSdfLabel = NULL;
                    ip->pSdfValue = NULL;
                    ip->nInputType = INPUT_MOLFILE;
                }
            }
            else if (INPUT_NONE == ip->nInputType && !inchi_stricmp(pArg, "MOL"))
            {
                ip->nInputType = INPUT_MOLFILE;
            }
            else if (INPUT_NONE == ip->nInputType && !inchi_stricmp(pArg, "SDF"))
            {
                ip->nInputType = INPUT_MOLFILE;
            }
            /*--- Output options ---*/
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
            else if (!inchi_stricmp(pArg, "Tabbed") || !inchi_stricmp(pArg, "Tab"))
            {
                bOutputStyle |= INCHI_OUT_TABBED_OUTPUT;
            }
#endif

            /* Removed condition
                #if ( defined(BUILD_WITH_ENG_OPTIONS) || defined(TARGET_LIB_FOR_WINCHI) )
                which is always true, as we already are under condition
                'if ( is_gui && ...)' and is_gui==1 means TARGET_LIB_FOR_WINCHI
            */
            /* #if ( defined(BUILD_WITH_ENG_OPTIONS) || defined(TARGET_LIB_FOR_WINCHI) )  */
            else if (!inchi_stricmp(pArg, "SDFID"))
            {
                ip->bGetSdfileId = 1;
            }
            else if (!inchi_stricmp(pArg, "PLAIN"))
            {
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT;
            }
            else if (!inchi_stricmp(pArg, "ANNPLAIN"))
            {
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT_COMMENTS;
            }
            else if (!inchi_memicmp(pArg, "AUXINFO:", 8) && isdigit(UCINT pArg[8]))
            {
                k = strtol(pArg + 8, NULL, 10);
                if (k == 0)
                {
                    bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
                    bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
                }
                else if (k == 1)
                {
                    bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include full aux info */
                }
                else if (k == 2)
                {
                    bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                    bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
                }
                else
                {
                    bINChIOutputOptions = k;  /* override everything */
                }
            }
            else if (!inchi_stricmp(pArg, "MERGE"))
            {
                bMergeAllInputStructures = 1;
            }

            else if (!inchi_stricmp(pArg, "PGO"))
            {
                ip->bSaveAllGoodStructsAsProblem = 1;
            }
            else if (!inchi_stricmp(pArg, "DCR"))
            {
                bDisplayCompositeResults = 1;
            }

            else if (!inchi_stricmp(pArg, "DSB"))
            {
                nMode |= REQ_MODE_NO_ALT_SBONDS;
            }
            else if (!inchi_stricmp(pArg, "NOHDR"))
            {
                bNoStructLabels = 1;
            }
            else if (!inchi_stricmp(pArg, "NoVarH"))
            {
                bTgFlagVariableProtons = 0;
            }
            /*--- (hidden) Old structure-perception and InChI creation options ---*/
            /*--- (engineering) Old structure-perception and InChI creation options ---*/
            else if (!inchi_stricmp(pArg, "NOUUSB"))
            {
                nMode |= REQ_MODE_SB_IGN_ALL_UU;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "NOUUSC"))
            {
                nMode |= REQ_MODE_SC_IGN_ALL_UU;
                bStdFormat = 0;
            }
#if ( FIX_ADJ_RAD == 1 )
            else if (!inchi_stricmp(pArg, "FixRad"))
            {
                bFixAdjacentRad = 1;
                bStdFormat = 0;
            }
#endif

#if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
            else if (!inchi_stricmp(pArg, "TestRenum") && developer_options)
            {
                ip->bRenumber = 1;
            }
#endif

#if ( UNDERIVATIZE == 1 )
            else if (!inchi_stricmp(pArg, "DoDRV"))
            {
                ip->bUnderivatize = 1;
                bStdFormat = 0;
            }
#if( UNDERIVATIZE_REPORT == 1 )
            else if (!inchi_stricmp(pArg, "DoDrvReport"))
            {
                ip->bUnderivatize = 3;
                bStdFormat = 0;
            }
#endif
#endif
#if ( RING2CHAIN == 1 )
            else if (!inchi_stricmp(pArg, "DoR2C"))
            {
                ip->bRing2Chain = 1;
                bStdFormat = 0;
            }
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
            else if (!inchi_stricmp(pArg, "DoneOnly"))
            {
                ip->bIgnoreUnchanged = 1;
                bStdFormat = 0;
            }
#endif
            else if (!inchi_stricmp(pArg, "NoADP"))
            {
                bTgFlagHardAddRenProtons = 0;
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "DISCONSALT:", 11))
            {
                bDisconnectSalts = (0 != strtol(pArg + 11, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "DISCONMETAL:", 12))
            {
                bDisconnectCoord = (0 != strtol(pArg + 12, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "RECONMETAL:", 11))
            {
                bReconnectCoord = (0 != strtol(pArg + 11, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "DISCONMETALCHKVAL:", 18))
            {
                bDisconnectCoordChkVal = (0 != strtol(pArg + 18, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "MOVEPOS:", 8))
            {
                bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "MERGESALTTG:", 12))
            {
                bMergeSaltTGroups = (0 != strtol(pArg + 12, NULL, 10));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "UNCHARGEDACIDS:", 15))
            {
                bUnchargedAcidTaut = (0 != strtol(pArg + 15, NULL, 16));
                bStdFormat = 0;
            }
            else if (!inchi_memicmp(pArg, "ACIDTAUT:", 9))
            {
                bAcidTautomerism = c = (int)strtol(pArg + 9, NULL, 10);
                if (0 <= c && c <= 2)  bAcidTautomerism = c;
                /*else bNotRecognized = 2*bReleaseVersion;*/
                bStdFormat = 0;
            }
            /*--- (hidden) Old output and other options ---*/
            else if (!inchi_memicmp(pArg, "O:", 2))
            {
                bNameSuffix = 1;
                strncpy(szNameSuffix, pArg + 2, sizeof(szNameSuffix) - 1);
            }
            else if (!inchi_memicmp(pArg, "OP:", 3))
            {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg + 3, sizeof(szOutputPath) - 1);
                }
            else if (!inchi_stricmp(pArg, "ALT"))
            {
                ip->bAbcNumbers = 1;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "SCT"))
            {
                ip->bCtPredecessors = 1;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "CMP"))
            {
                bCompareComponents = CMP_COMPONENTS;
            }
            else if (!inchi_stricmp(pArg, "CMPNONISO"))
            {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            }
            else if (!inchi_stricmp(pArg, "PW"))
            {
                ip->bSaveWarningStructsAsProblem = 1;
            }

            else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
            {
                mdbr = (int)strtol(pArg + 4, NULL, 10);
            }
            else if (!inchi_memicmp(pArg, "DISCONSALT:", 11) && developer_options)
            {
                bDisconnectSalts = (0 != strtol(pArg + 11, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "DISCONMETAL:", 12) && developer_options)
            {
                bDisconnectCoord = (0 != strtol(pArg + 12, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "RECONMETAL:", 11) && developer_options)
            {
                bReconnectCoord = (0 != strtol(pArg + 11, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "DISCONMETALCHKVAL:", 18) && developer_options)
            {
                bDisconnectCoordChkVal = (0 != strtol(pArg + 18, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
            {
                bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "MERGESALTTG:", 12) && developer_options)
            {
                bMergeSaltTGroups = (0 != strtol(pArg + 12, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "UNCHARGEDACIDS:", 15) && developer_options)
            {
                bUnchargedAcidTaut = (0 != strtol(pArg + 15, NULL, 16));;
            }
            else if (!inchi_memicmp(pArg, "ACIDTAUT:", 9) && developer_options)
            {
                bAcidTautomerism = c = (int)strtol(pArg + 9, NULL, 10);
                if (0 <= c && c <= 2)
                {
                    bAcidTautomerism = c;
                }
                /*else bNotRecognized = 2*bReleaseVersion;*/
            }

            else
            {
                /*for ( k = 0; c=pArg[k]; k ++ )*/
                k = 0;
                c = pArg[k]; /* prohibit multiple option concatenations, strict syntax check 2008-11-05 DT  */
                {
                    c = toupper(c);
                    switch (c)
                    {
                    case 'D':
                        bDisplay |= 1;
                        if ((pArg[k + 1] == 'C' || pArg[k + 1] == 'c') && !pArg[k + 2])
                        {
                            bDisplay |= 1;
                            k++;
                            ip->bDisplayEachComponentINChI = 1;
                        }
                        else if (!pArg[k + 1])
                        {
                            bDisplay |= 1;
                        }
                        break;
                    case 'W':
                        if (pArg[k + 1] == 'D')
                        {
                            /* restore Display Time functionality */
                            c = 'D';
                            k++;
                        }
                        t = strtod(pArg + k + 1, (char**)&q); /*  cast deliberately discards 'const' qualifier */
                        if ((q > pArg + k + 1 && errno == ERANGE) || t < 0.0 || t * 1000.0 >(double)ULONG_MAX) /* djb-rwth: addressing LLVM warning */
                        {
                            ul = 0;
                        }
                        else
                        {
                            ul = (unsigned long)(t * 1000.0);
                        }
                        if ( /*q > pArg+k &&*/ !*q)
                        {
                            k = q - pArg - 1; /* k will be incremented by the for() cycle */
                            switch (c)
                            {
                            case 'D':
                                *ulDisplTime = ul;
                                break;
                            case 'W':
                                ip->msec_MaxTime = ul;
                                break;
                            }
                        }
                        break;
                    case 'F':
                        c = (int)strtol(pArg + k + 1, (char**)&q, 10); /*  cast deliberately discards 'const' qualifier */
                        if (q > pArg + k && !*q)
                        {
                            k = q - pArg - 1;
                            if (abs(c) > 5)
                            {
                                nFontSize = -c;  /* font size 5 or less is too small */
                            }
                        }
                        break;
                    default:
                        if (!pArg[k + 1])
                        {
                            switch (c)
                            {
                            case 'B':
                                nMode |= REQ_MODE_BASIC;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                bStdFormat = 0;
                                break;
                            case 'T':
                                nMode |= REQ_MODE_TAUT;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'I':
                                nMode |= REQ_MODE_ISO;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'N':
                                nMode |= REQ_MODE_NON_ISO;
                                nReleaseMode = 0;
                                bStdFormat = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'S':
                                nMode |= REQ_MODE_STEREO;
                                nReleaseMode = 0;
                                /*bNotRecognized = bReleaseVersion;*/
                                break;
                            case 'E':
                                if (nReleaseMode & REQ_MODE_STEREO)
                                {
                                    nReleaseMode ^= REQ_MODE_STEREO;
                                    bStdFormat = 0;
                                }
                                break;

#ifndef TARGET_LIB_FOR_WINCHI
                            default:
                                inchi_ios_eprint(log_file, "Unrecognized optionQ1: \"%c\".\n", c);

#endif
                            }
                        }


#ifndef TARGET_LIB_FOR_WINCHI
                        else
                        {
                            inchi_ios_eprint(log_file, "Unrecognized optionQ2: \"%c\".\n", c);
                        }
#endif
                    }
                    /*
                        if ( bNotRecognized && bNotRecognized == bReleaseVersion ) {
                        inchi_ios_eprint(stderr, "Unrecognized option: \"%c\".\n", c);
                        bNotRecognized = 0;
                        }
                    */
                }
            }

            /*
                if ( bNotRecognized && bNotRecognized == 2*bReleaseVersion )
                {
                inchi_ios_eprint(stderr, "Unrecognized option: \"%s\".\n", argv[i]);
                bNotRecognized = 0;
                }
            */

        } /* eof Parsing TARGET_LIB_FOR_WINCHI GUI options (and v. 0.9xx Beta as well)  */

        else if ((bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1])
        {
            /* Parsing stand-alone executable/libinchi options */

            pArg = argv[i] + 1;

            bRecognizedOption = 2;
            bVer1Options += 2;
            /* always on: REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO */

            got = set_common_options_by_parg(pArg, developer_options, ip, &bVer1DefaultMode, &nMode,
                &bINChIOutputOptions, &bINChIOutputOptions2,
                &bStdFormat, &bHashKey, &bHashXtra1, &bHashXtra2,
                &bFixSp3bug, &bFixFB2,
                &bAddPhosphineStereo, &bAddArsineStereo,
                &bNoStructLabels, &bPointedEdgeStereo,
                &bDoNotAddH, &bForcedChiralFlag, &bReconnectCoord,
                &bKetoEnolTaut,
                &b15TautNonRing,
                &bPT_06_00_Taut, &bPT_13_00_Taut, &bPT_16_00_Taut,
                &bPT_18_00_Taut, &bPT_22_00_Taut, &bPT_39_00_Taut,
                &bLooseTSACheck,
                &bLargeMolecules, &bPolymers,
                &bFoldPolymerSRU, &bFrameShiftScheme,
                &bStereoAtZz, &bNPZz,
                &bNoWarnings, &bMergeHash, &bHideInChI);

            if (got)
            {
                ;
            }
            /* Input */
            else if (!inchi_stricmp(pArg, "STDIO"))
            {
                bNameSuffix = 0;
            }
            else if ( /* INPUT_NONE == ip->nInputType &&*/
                !inchi_memicmp(pArg, "SDF:", 4))
            {
                /* SDfile label */
                k = 0;
                mystrncpy(ip->szSdfDataHeader, pArg + 4, MAX_SDF_HEADER + 1);
                lrtrim(ip->szSdfDataHeader, &k);
                if (k)
                {
                    ip->pSdfLabel = ip->szSdfDataHeader;
                    ip->pSdfValue = szSdfDataValue;
                    if (INPUT_NONE == ip->nInputType)
                    {
                        ip->nInputType = INPUT_SDFILE;
                    }
                }
                else
                {
                    ip->pSdfLabel = NULL;
                    ip->pSdfValue = NULL;
                    if (INPUT_NONE == ip->nInputType)
                    {
                        ip->nInputType = INPUT_MOLFILE;
                    }
                }
            }

            else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
            {
                mdbr = (int)strtol(pArg + 4, NULL, 10);
            }

            /* Output */
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
            else if (!inchi_stricmp(pArg, "Tabbed") || !inchi_stricmp(pArg, "Tab"))
            {
                bOutputStyle |= INCHI_OUT_TABBED_OUTPUT;
            }
#endif
            else if (!inchi_stricmp(pArg, "OUTPUTSDF"))
            {
                /* output SDfile */
                bOutputMolfileOnly = 1;
            }
            else if (!inchi_stricmp(pArg, "SdfAtomsDT"))
            {
                /* output isotopes H as D and T in SDfile */
                bOutputMolfileDT = 1;
            }
            else if (!inchi_stricmp(pArg, "D"))
            {
                /* display the structures */
                bDisplay |= 1;
            }
            else if (!inchi_memicmp(pArg, "F", 1) && (c = (int)strtol(pArg + 1, (char**)&q, 10), q > pArg + 1))
            {
                nFontSize = -c;                      /* struct. display font size */
            }
            else if (!inchi_stricmp(pArg, "EQU"))
            {
                bCompareComponents = CMP_COMPONENTS;
            }
#if( OUTPUT_FILE_EXT == 1 )
            else if (pArg[0] == '.' && numOutNameExt < (int)(sizeof(szOutNameExt) / sizeof(szOutNameExt[0])))
            {
                strncpy(szOutNameExt[numOutNameExt], pArg, sizeof(szOutNameExt[0]) - 1);
                numOutNameExt++;
                if (ip->path[numOutNameExt])
                {
                    ip->path[numOutNameExt] = ""; /*strcpy( ip->path[numOutNameExt], "");*/
                }
            }
#endif
            /* djb-rwth: avoiding Error 98 for empty .mol files -- GH issue #25, thanks to @wijnand1 */
            else if (!inchi_stricmp(pArg, "WarnOnEmptyStructure"))
            {
                ip->bAllowEmptyStructure = 1;
            }

            /* Generation options */
            else if (!inchi_memicmp(pArg, "W", 1))
            {
                long timeout_value;
                const char c1 = *(pArg + 1);
                if (c1 && (c1 == 'M' || c1 == 'm'))
                {
                    /* "WMnumber", milliseconds */
                    timeout_value = strtol(pArg + 2, (char**)&q, 10);
                    if (timeout_value && q > pArg + 2 && *q == '\0')
                    {
                        if (errno == ERANGE || timeout_value < 0.0 || timeout_value>LONG_MAX) /* djb-rwth: addressing coverity ID #499550 -- the condition takes into account all possible overflows/errors */
                        {
                            timeout_value = 0;
                            timeout_set_warning = 1;
                        }
                        timeout_set_error = 0;
                    }
                    else
                    {
                        timeout_set_error = 1;
                    }
                }
                else
                {
                    /* expect "Wnumber", seconds */
                    t = strtod(pArg + 1, (char**)&q);
                    if (t && q > pArg + 1)
                    {
                        if (*q != '\0')
                        {
                            timeout_set_warning = 1;
                        }
                        if (errno == ERANGE || t < 0.0 || t * 1000.0 >(double)LONG_MAX)
                        {
                            timeout_value = 0;
                            timeout_set_warning = 1;
                        }
                        else
                        {
                            timeout_value = (long)(t * 1000.0);  /* max. time per structure */
                        }
                        timeout_set_error = 0;
                    }
                    else
                    {
                        timeout_set_error = 1;
                    }
                }
                if (timeout_set_error == 0)
                {
                    ip->msec_MaxTime = timeout_value;
                }
            }

            /*--- Conversion modes ---*/
#if ( READ_INCHI_STRING == 1 )

/*#if (BUILD_WITH_ENG_OPTIONS==1)*/
            else if (!inchi_stricmp(pArg, "InChI2InChI"))
            {
                /* Read InChI Identifiers and output InChI Identifiers */
                ip->nInputType = INPUT_INCHI;
                ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
                ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
            }
            /*#endif*/

            else if (!inchi_stricmp(pArg, "InChI2Struct"))
            {
                /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_TO_STRUCTURE;
                ip->bReadInChIOptions &= ~READ_INCHI_OUTPUT_INCHI;
                ip->nInputType = INPUT_INCHI;
            }

            else if (!inchi_stricmp(pArg, "KeepBalanceP") && developer_options)
            {
                /* When spliting InChI Identifiers into components: */
                /* If MobileH output then add p to each component;  */
                /* Otherwise add one more component containing balance */
                /* of protons and exchangeable isotopic H */
                ip->bReadInChIOptions |= READ_INCHI_KEEP_BALANCE_P;
                bStdFormat = 0;
            }
#endif

            /*--- (engineering) Undo bug/draw fixes options ---*/

            /* (developer_options) Old structure-perception and InChI creation options */
#if ( FIX_ADJ_RAD == 1 )
            else if (!inchi_stricmp(pArg, "FixRad") && developer_options)
            {
                bFixAdjacentRad = 1;
                bStdFormat = 0;
            }
#endif


#if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
            else if (!inchi_stricmp(pArg, "TestRenum") && developer_options)
            {
                ip->bRenumber = 1;
            }
#endif

#if ( UNDERIVATIZE == 1 )
            else if (!inchi_stricmp(pArg, "DoDRV") && developer_options)
            {
                ip->bUnderivatize = 1;
                bStdFormat = 0;
            }
#if( UNDERIVATIZE_REPORT == 1 )
            else if (!inchi_stricmp(pArg, "DoDrvReport"))
            {
                ip->bUnderivatize = 3;
                bStdFormat = 0;
            }
#endif
#endif
#if ( RING2CHAIN == 1 )
            else if (!inchi_stricmp(pArg, "DoR2C") && developer_options)
            {
                ip->bRing2Chain = 1;
                bStdFormat = 0;
            }
#endif
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
            else if (!inchi_stricmp(pArg, "DoneOnly") && developer_options)
            {
                ip->bIgnoreUnchanged = 1;
                bStdFormat = 0;
            }
#endif
            else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
            {
                bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
                bStdFormat = 0;
            }

            else if (!inchi_stricmp(pArg, "NoADP") && developer_options)
            {
                bTgFlagHardAddRenProtons = 0;
                bStdFormat = 0;
            }
            /* Tautomer perception off */
            else if (!inchi_stricmp(pArg, "EXACT") && developer_options)
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "ONLYRECSALT") && developer_options)
            {
                /* do not disconnect salts */
                bDisconnectSalts = 0;
                bStdFormat = 0;
            }
            else if ((!inchi_stricmp(pArg, "ONLYEXACT") || !inchi_stricmp(pArg, "ONLYFIXEDH")) && developer_options)
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;
                bVer1DefaultMode &= ~REQ_MODE_TAUT;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "ONLYNONISO") && developer_options)
            {
                bVer1DefaultMode |= REQ_MODE_NON_ISO;
                bVer1DefaultMode &= ~REQ_MODE_ISO;
                bStdFormat = 0;
            }
            else if (!inchi_stricmp(pArg, "TAUT") && developer_options)
            {
                bVer1DefaultMode &= ~REQ_MODE_BASIC;
                bVer1DefaultMode |= REQ_MODE_TAUT;
            }
            else if (!inchi_stricmp(pArg, "ONLYRECMET") && developer_options)
            {
                /* do not disconnect metals */
                bDisconnectCoord = 0;
                bStdFormat = 0;
            }

            /*--- (hidden) Old output and other options ---*/

            else if (!inchi_stricmp(pArg, "SdfSplit") && developer_options)
            {
                /* Split single Molfiles into disconnected components */
                bOutputMolfileSplit = 1;
            }
            else if (!inchi_stricmp(pArg, "DCR") && developer_options)
            {
                bDisplayCompositeResults = 1;
            }
            else if ((!inchi_stricmp(pArg, "AUXFULL") || !inchi_stricmp(pArg, "AUXMAX")) && developer_options)
            {
                /* full aux info */
                bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include short aux info */
            }
            else if (!inchi_stricmp(pArg, "AUXMIN") && developer_options)
            {
                /* minimal aux info */
                bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
                bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
            }

#if ( READ_INCHI_STRING == 1 )
            else if (!inchi_stricmp(pArg, "DDSRC") && developer_options)
            {
                bDisplayIfRestoreWarnings = 1;  /* InChI->Structure debugging: Display Debug Structure Restore Components */
            }
#endif

            else if (!inchi_stricmp(pArg, "NoVarH") && developer_options)
            {
                bTgFlagVariableProtons = 0;
            }
            else if (!inchi_stricmp(pArg, "FULL") && developer_options)
            {
                bVer1DefaultMode = VER103_DEFAULT_MODE;
                nMode = 0;
                bReconnectCoord = 1;            /* full output */
                bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0) | INCHI_OUT_SHORT_AUX_INFO;
                ip->bCtPredecessors = 0;
                ip->bAbcNumbers = 0;
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
            }
            else if (!inchi_stricmp(pArg, "MIN") && developer_options)
            {
                bVer1DefaultMode = VER103_DEFAULT_MODE;
                nMode = 0;
                bReconnectCoord = 1;            /* minimal output */
                bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0) | INCHI_OUT_NO_AUX_INFO;            /* minimal compressed output */
                ip->bCtPredecessors = 1;
                ip->bAbcNumbers = 1;
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
            }
            else if (!inchi_stricmp(pArg, "COMPRESS") && developer_options)
            {
                ip->bAbcNumbers = 1;
                ip->bCtPredecessors = 1;             /* compressed output */
            }

#if ( READ_INCHI_STRING == 1 )
            else if (!inchi_stricmp(pArg, "InChI2InChI")) /*&& developer_options)*/
            {
                /* Read InChI Identifiers and output InChI Identifiers */
                ip->nInputType = INPUT_INCHI;
                ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
                ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
            }

            else if (!inchi_stricmp(pArg, "SplitInChI") && developer_options)
            {
                /* Split InChI Identifiers into components */
                ip->bReadInChIOptions |= READ_INCHI_SPLIT_OUTPUT;
            }
#endif

            else if (!inchi_stricmp(pArg, "MOLFILENUMBER") && developer_options)
            {
                ip->bGetMolfileNumber |= 1;
            }
            else if (!inchi_stricmp(pArg, "OutputPLAIN") && developer_options)
            {
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT;
            }
            else if (!inchi_stricmp(pArg, "OutputANNPLAIN") && developer_options)
            {
                bOutputStyle |= INCHI_OUT_PLAIN_TEXT_COMMENTS;
                bOutputStyle |= INCHI_OUT_WINCHI_WINDOW; /* debug */
            }
            else if ((!inchi_stricmp(pArg, "ONLYEXACT") || !inchi_stricmp(pArg, "ONLYFIXEDH")) && developer_options)
            {
                bVer1DefaultMode |= REQ_MODE_BASIC;
                bVer1DefaultMode &= ~REQ_MODE_TAUT;
            }
            else if (!inchi_stricmp(pArg, "ONLYNONISO") && developer_options)
            {
                bVer1DefaultMode |= REQ_MODE_NON_ISO;
                bVer1DefaultMode &= ~REQ_MODE_ISO;
            }
            else if (!inchi_stricmp(pArg, "TAUT") && developer_options)
            {
                bVer1DefaultMode &= ~REQ_MODE_BASIC;
                bVer1DefaultMode |= REQ_MODE_TAUT;
            }
            else if (!inchi_stricmp(pArg, "ONLYRECMET") && developer_options)
            {  /* do not disconnect metals */
                bDisconnectCoord = 0;
            }
            else if (!inchi_stricmp(pArg, "ONLYRECSALT") && developer_options)
            {  /* do not disconnect salts */
                bDisconnectSalts = 0;
            }
            else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
            {   /* added -- 2010-03-01 DT */
                bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
            }
            else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
            {
                mdbr = (int)strtol(pArg + 4, NULL, 10);
            }
            else if (!inchi_stricmp(pArg, "EQU") && developer_options)
            {
                bCompareComponents = CMP_COMPONENTS;
            }
            else if (!inchi_stricmp(pArg, "EQUNONISO") && developer_options)
            {
                bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
            }
            else if (!inchi_memicmp(pArg, "OP:", 3) && developer_options)
            {
                bOutputPath = 1;
                strncpy(szOutputPath, pArg + 3, sizeof(szOutputPath) - 1);
            }
            /* eof developer_options */

            /* Display unrecognized option */
            else
            {
                bRecognizedOption = 0;

#ifndef TARGET_LIB_FOR_WINCHI
                inchi_ios_eprint(log_file, "Unrecognized optionQ3: \"%s\".\n", pArg);
#endif
            }
            bVer1Options |= bRecognizedOption;

        } /* eof Parsing stand-alone executable/libinchi options */

        else if (ip->num_paths < MAX_NUM_PATHS)
        {
            char* sz;
#if( ALLOW_EMPTY_PATHS == 1 )
            if (argv[i])
#else
            if (argv[i] && argv[i][0])
#endif
            {
                if ((sz = (char*)inchi_malloc((strlen(argv[i]) + 1) * sizeof(sz[0])))) /* djb-rwth: addressing LLVM warning */
                {
                    strcpy(sz, argv[i]);
                }

                ip->path[ip->num_paths++] = sz;
            }
        }
    } /*  eof parsing argv loop */


    /* Print messages and set controil variables according to just parsed options */

    /* Timeout option(s) */
    if (timeout_set_warning)
    {
        inchi_ios_eprint(log_file, "Warning: timeout value may have been modified (truncated?) due to number formatting issues;\n");
    }
    if (timeout_set_error)
    {
        inchi_ios_eprint(log_file, "Warning: specified timeout value was ignored due to invalid number format, using the default;\n");
    }

    /* InChIKey option(s) */
    if (bHashKey != 0)
    {
        /* Suppress InChIKey calculation if:                */
        /* compressed output OR Inchi2Struct OR Inchi2Inchi */
        if ((ip->bAbcNumbers == 1) && (ip->bCtPredecessors == 1))
        {
            bHashKey = 0;
#ifndef TARGET_LIB_FOR_WINCHI
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'Compress' option\n");
            return -1;
#endif
        }
        if (ip->nInputType == INPUT_INCHI)
        {
            bHashKey = 0;
#ifndef TARGET_LIB_FOR_WINCHI
            inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available in InChI conversion mode\n");
            return -1;
#endif
        }
        else
        {
            if (bOutputMolfileOnly == 1)
            {
                bHashKey = 0;
#ifndef TARGET_LIB_FOR_WINCHI
                inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n");
                return -1;
#endif
            }
        }
    }

    if (bNameSuffix || bOutputPath)
    {
        const char szNUL[] = "NUL";
        /* fix for AMD processor: use const char[] instead of just "NUL" constant 2008-11-5 DT */
        const char* p = NULL;
        char* p_prev = NULL; /* djb-rwth: avoiding use of memory after it is freed */
        char* r = NULL;
        char* sz;
        int len;
        /*  find the 1st path */
        for (i = 0; i < MAX_NUM_PATHS; i++)
        {
            if (!p && ip->path[i] && ip->path[i][0])
            {
                p = ip->path[i];
                break;
            }
        }
        /* fix output path */
        if (bOutputPath && szOutputPath[0] && p)
        {
            /* remove last slash */
            len = (int)strlen(szOutputPath);
            if (len > 0 && szOutputPath[len - 1] != INCHI_PATH_DELIM)
            {
                szOutputPath[len++] = INCHI_PATH_DELIM;
                szOutputPath[len] = '\0';
            }
            if (len > 0 && (r = (char*)strrchr(p, INCHI_PATH_DELIM)) && r[1])
            {
                strcat(szOutputPath, r + 1);
                p = szOutputPath;
            }
        }
/* djb-rwth: copying the value of p */
#ifdef _WIN32
        p_prev = _strdup(p);
#else
        p_prev = inchi__strdup(p);
#endif
        /*  add missing paths */
        /* djb-rwth: this whole block had to be rewritten to avoid use of memory after it is freed */
        for (i = 0; p_prev && i < MAX_NUM_PATHS; i++)
        {
            /* fix for AMD processor: changed order 2008-11-5 DT */
            if (!ip->path[i] || !ip->path[i][0])
            {
#if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
                char* pLastExt = (i && numOutNameExt >= i) ? strrchr((char*)p_prev, '.') : 0;
                char* pLastSlash = (i && numOutNameExt >= i) ? strrchr((char*)p_prev, INCHI_PATH_DELIM) : 0;
                if (pLastExt && pLastSlash && pLastSlash > pLastExt)
                    pLastExt = NULL;
#else
                char* pLastExt = NULL;
#endif
                len = (int)strlen(p_prev) + strlen(szNameSuffix) + strlen(ext[i]);
                if ((sz = (char*)inchi_malloc(((long long)len + 1) * sizeof(sz[0])))) /* djb-rwth: cast operators added; addressing and ignoring LLVM warnings */
                {
                    strcpy(sz, p_prev); /* djb-rwth: fix for use of memory after being freed */
#if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
                    if (pLastExt)
                    {
                        strcpy(sz + (pLastExt - p), szOutNameExt[i - 1]);
                    }
#endif
                    strcat(sz, szNameSuffix);
                    if (!pLastExt)
                        strcat(sz, ext[i]);
                    ip->num_paths++;
                    if (ip->path[i])
                    {
                        inchi_free((char*)ip->path[i]); /* eliminate memory leak 2013-12-18 DCh */
                    }
                    ip->path[i] = sz;
                }
            }
            else if (!inchi_stricmp(ip->path[i], szNUL))
            {
                inchi_free((char*)ip->path[i]); /* cast deliberately const qualifier */
                ip->path[i] = NULL;
            }
        }
        free(p_prev); /* djb-rwth: freeing memory reserved for auxiliary variable */
    }

    /* inchi2inchi and inchi2struct option(s) */
#if ( READ_INCHI_STRING == 1 )
    if (INPUT_INCHI == ip->nInputType)
    {
        bCompareComponents = 0;
        /*bDisplayCompositeResults = 0;*/
#if ( I2S_MODIFY_OUTPUT == 1 )
        if (!(ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE))
#endif
        {
            bOutputMolfileOnly = 0;
            bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;
            bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
            bINChIOutputOptions &= ~INCHI_OUT_ONLY_AUX_INFO;
            /* bNoStructLabels   = 1; */
        }
        ip->bDisplayIfRestoreWarnings = bDisplayIfRestoreWarnings;
        if (!(bINChIOutputOptions &
            (INCHI_OUT_SDFILE_ONLY |       /* not in bINChIOutputOptions yet */
                INCHI_OUT_PLAIN_TEXT |        /* not in bINChIOutputOptions yet */
                INCHI_OUT_PLAIN_TEXT_COMMENTS /* not in bINChIOutputOptions yet */
                )
            )
#if ( I2S_MODIFY_OUTPUT == 1 )
            && !bOutputMolfileOnly
            && !(bOutputStyle & (INCHI_OUT_PLAIN_TEXT))
#endif
            )
        {
            bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT;
        }
    }
#endif

    if (bVer1Options)
    {
        nMode |= bVer1DefaultMode;
    }
    else if (bReleaseVersion)
    {
        nMode |= nReleaseMode;
    }

#if ( defined(COMPILE_ANSI_ONLY) || defined(TARGET_LIB_FOR_WINCHI) )
    if (bCompareComponents && !(bDisplay & 1))
    {
        bCompareComponents = 0;
    }
#endif

    /* Save original options */
    /* nOrigMode = nMode; */
#ifndef COMPILE_ANSI_ONLY
    ip->dp.sdp.nFontSize = nFontSize;
    ip->dp.sdp.ulDisplTime = *ulDisplTime;
    ip->bDisplay = bDisplay;
#ifdef TARGET_LIB_FOR_WINCHI
    ip->bDisplayCompositeResults = bDisplay;
#else
    ip->bDisplayCompositeResults = bDisplayCompositeResults;
#endif
#else
    ip->bDisplayEachComponentINChI = 0;
    bCompareComponents = 0;
#endif

    ip->bMergeAllInputStructures = bMergeAllInputStructures;
    ip->bDoNotAddH = bDoNotAddH;

    /*  Set default options */
    if (!nMode || nMode == REQ_MODE_STEREO)
    {
        /*  requested all output */
        nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO);
    }
    else
    {
        if (!(nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)))
        {
            nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT);
        }
        if ((nMode & REQ_MODE_STEREO) && !(nMode & (REQ_MODE_ISO | REQ_MODE_NON_ISO)))
        {
            nMode |= (REQ_MODE_ISO | REQ_MODE_NON_ISO);
        }
    }
    /*  if the user requested isotopic then unconditionally add non-isotopic output. */
    if (nMode & REQ_MODE_ISO)
    {
        nMode |= REQ_MODE_NON_ISO;
    }

#if ( MIN_SB_RING_SIZE > 0 )
    if (mdbr)
    {
        nMinDbRinSize = mdbr;
    }
    nMode |= (nMinDbRinSize << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK;
#endif

    /* Input file */
    if (ip->nInputType == INPUT_NONE && ip->num_paths > 0)
    {
        ip->nInputType = INPUT_MOLFILE; /*  default */
    }
    ip->nMode = nMode;
    /* Compare components */
    if ((bCompareComponents & CMP_COMPONENTS) && (nMode & REQ_MODE_BASIC))
    {
        bCompareComponents |= CMP_COMPONENTS_NONTAUT; /* compare non-tautomeric */
    }
    ip->bCompareComponents = bCompareComponents;
    /* Output */
    ip->bINChIOutputOptions = bINChIOutputOptions | (bOutputMolfileOnly ? INCHI_OUT_SDFILE_ONLY : 0);
    if (bOutputMolfileOnly)
    {
        bOutputStyle &= ~(INCHI_OUT_PLAIN_TEXT |
            INCHI_OUT_PLAIN_TEXT_COMMENTS |
            INCHI_OUT_TABBED_OUTPUT);
#if ( SDF_OUTPUT_DT == 1 )
        ip->bINChIOutputOptions |= bOutputMolfileDT ? INCHI_OUT_SDFILE_ATOMS_DT : 0;
        ip->bINChIOutputOptions |= bOutputMolfileSplit ? INCHI_OUT_SDFILE_SPLIT : 0;
#endif
    }
#ifdef TARGET_LIB_FOR_WINCHI
    if (!(bDisplay & 1))
    {
        bOutputStyle &= ~(INCHI_OUT_PLAIN_TEXT_COMMENTS); /* do not ouput comments in wINChI text file results */
    }
    else
    {
        bOutputStyle |= INCHI_OUT_WINCHI_WINDOW;
    }
#endif
    ip->bINChIOutputOptions |= bOutputStyle;
    ip->bNoStructLabels = bNoStructLabels;

    /* Processing options */
    if (bForcedChiralFlag)
    {
        ip->bChiralFlag = bForcedChiralFlag;
    }

    /* Tautomeric/salts options */
    ip->bTautFlags = 0;
    ip->bTautFlagsDone = 0;

    /* Find regular tautomerism */
    ip->bTautFlags |= TG_FLAG_TEST_TAUT__ATOMS;

    /* Disconnect salts */
    ip->bTautFlags |= bDisconnectSalts ? TG_FLAG_DISCONNECT_SALTS : 0;

    /* If possible, find long-range H/(-) taut. on =C-OH, >C=O    */
    ip->bTautFlags |= bAcidTautomerism ? TG_FLAG_TEST_TAUT__SALTS : 0;

    /* Allow long-range movement of N(+), P(+) charges           */
    ip->bTautFlags |= bMovePositiveCharges ? TG_FLAG_MOVE_POS_CHARGES : 0;

    /* Multi-attachement long-range H/(-) taut. on =C-OH, >C=O   */
    ip->bTautFlags |= (bAcidTautomerism > 1) ? TG_FLAG_TEST_TAUT2_SALTS : 0;

    /* (Debug) allow to find long-range H-only tautomerism on =C-OH, >C=O */
    ip->bTautFlags |= (bUnchargedAcidTaut == 1) ? TG_FLAG_ALLOW_NO_NEGTV_O : 0;

    /* Merge =C-OH and >C=O containing t-groups and other =C-OH groups */
    ip->bTautFlags |= bMergeSaltTGroups ? TG_FLAG_MERGE_TAUT_SALTS : 0;
    ip->bTautFlags |= bDisconnectCoord ? TG_FLAG_DISCONNECT_COORD : 0;
    ip->bTautFlags |= (bDisconnectCoord && bReconnectCoord) ? TG_FLAG_RECONNECT_COORD : 0;
    ip->bTautFlags |= bDisconnectCoordChkVal ? TG_FLAG_CHECK_VALENCE_COORD : 0;
    ip->bTautFlags |= bTgFlagVariableProtons ? TG_FLAG_VARIABLE_PROTONS : 0;
    ip->bTautFlags |= bTgFlagHardAddRenProtons ? TG_FLAG_HARD_ADD_REM_PROTONS : 0;
    ip->bTautFlags |= bKetoEnolTaut ? TG_FLAG_KETO_ENOL_TAUT : 0;
    ip->bTautFlags |= b15TautNonRing ? TG_FLAG_1_5_TAUT : 0;

    /*^^^ IPl 2020-04-02 added to forcefully enable new tauto rules (in test_ixa, etc.) */
    /*^^^ IPl 2020-10-26 removed set to 1
    bPT_22_00_Taut = 1;
    bPT_16_00_Taut = 1;
    bPT_06_00_Taut = 1;
    bPT_39_00_Taut = 1;
    */

    /*^^^ IPl 2020-04-02 */
    ip->bTautFlags |= bPT_22_00_Taut ? TG_FLAG_PT_22_00 : 0;
    ip->bTautFlags |= bPT_16_00_Taut ? TG_FLAG_PT_16_00 : 0;
    ip->bTautFlags |= bPT_06_00_Taut ? TG_FLAG_PT_06_00 : 0;
    ip->bTautFlags |= bPT_39_00_Taut ? TG_FLAG_PT_39_00 : 0;
    ip->bTautFlags |= bPT_13_00_Taut ? TG_FLAG_PT_13_00 : 0;
    ip->bTautFlags |= bPT_18_00_Taut ? TG_FLAG_PT_18_00 : 0;

#ifdef STEREO_WEDGE_ONLY
    ip->bTautFlags |= bPointedEdgeStereo ? TG_FLAG_POINTED_EDGE_STEREO : 0;
#endif
#if ( FIX_ADJ_RAD == 1 )
    ip->bTautFlags |= bFixAdjacentRad ? TG_FLAG_FIX_ADJ_RADICALS : 0;
#endif
    ip->bTautFlags |= bAddPhosphineStereo ? TG_FLAG_PHOSPHINE_STEREO : 0;
    ip->bTautFlags |= bAddArsineStereo ? TG_FLAG_ARSINE_STEREO : 0;
    ip->bTautFlags |= bFixSp3bug ? TG_FLAG_FIX_SP3_BUG : 0;

    /* Bug fixes */
    if (bFixFB2)
    {
#if ( FIX_ISO_FIXEDH_BUG == 1 )
        ip->bTautFlags |= TG_FLAG_FIX_ISO_FIXEDH_BUG; /* accomodate FIX_ISO_FIXEDH_BUG */
#endif

#if ( FIX_TERM_H_CHRG_BUG == 1 )
        ip->bTautFlags |= TG_FLAG_FIX_TERM_H_CHRG_BUG; /* accomodate FIX_TERM_H_CHRG_BUG */
#endif

#if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 )
        ip->bINChIOutputOptions |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG;
#endif
    }

    if (!ip->nInputType)
    {
        ip->nInputType = INPUT_MOLFILE;
    }

    /* Check if /SNon requested turn OFF SUU/SLUUD */
    if (!(ip->nMode & REQ_MODE_STEREO))
    {
        ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        ip->nMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
    }

    /* Standard InChI ? */
    if (bStdFormat)
    {
        ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
    }

    /* InChIKey ? */
    if (!bHashKey)
    {
        ip->bCalcInChIHash = INCHIHASH_NONE;
    }
    else
    {
        ip->bCalcInChIHash = INCHIHASH_KEY;
    }
    /* Extension(s) to hash (in non-std mode only) ? */
    if (!bHashKey)
    {
        if ((bHashXtra1 != 0) || (bHashXtra2 != 0))
        {
            inchi_ios_eprint(log_file, "Hash extension(s) not generated: InChIKey not requested");
        }
    }
    else
    {
        if (bHashXtra1)
        {
            if (bHashXtra2)
            {
                ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1_XTRA2;
            }
            else
            {
                ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1;
            }
        }
        else if (bHashXtra2)
        {
            ip->bCalcInChIHash = INCHIHASH_KEY_XTRA2;
        }
    }

    ip->bLargeMolecules = bLargeMolecules;
    ip->bLooseTSACheck = bLooseTSACheck;

    ip->bNPZz = bNPZz;
    ip->bStereoAtZz = bStereoAtZz;  /*STEREO_AT_ZZ;*/

    ip->bNoWarnings = bNoWarnings;
    ip->bMergeHash = bMergeHash;
    ip->bHideInChI = bHideInChI;

    ip->bPolymers = bPolymers;
    ip->bFoldPolymerSRU = bFoldPolymerSRU;
    ip->bFrameShiftScheme = bFrameShiftScheme;

#ifdef TARGET_LIB_FOR_WINCHI
    ip->bLargeMolecules = 1;
    /*ip->bPolymers = POLYMERS_MODERN; */
    /*ip->bNPZz = 1;*/
    /*ip->bFrameShiftScheme = FSS_NONE;*/
    /*ip->bFoldPolymerSRU = 1;*/

#if ( UNDERIVATIZE == 1 )
    ip->bUnderivatize = 1;
    if (ip->bUnderivatize)
    {
        ip->bINChIOutputOptions &= ~INCHI_OUT_STDINCHI;
    }
#endif
    ip->bNoWarnings = 0;
    ip->bMergeHash = 0;
    ip->bHideInChI = 0;

#endif

    ip->bINChIOutputOptions2 = bINChIOutputOptions2;

    return 0;
}


/****************************************************************************
 Print input parameters
****************************************************************************/
int PrintInputParms(INCHI_IOSTREAM* log_file,
    INPUT_PARMS* ip)
{

#if (BUILD_WITH_ENG_OPTIONS==1)
    const int developer_options = 1;
#else
    const int developer_options = 0;
#endif
#ifdef TARGET_LIB_FOR_WINCHI
    int bInChI2Struct = 0; /* winchi-1 can not convert InChI to structure */
#else
    int bInChI2Struct = (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE) && ip->nInputType == INPUT_INCHI;
#endif

    INCHI_MODE nMode = ip->nMode;
    int k;
    int bStdFormat = 1;
    int first = 1;


    if (!(ip->bINChIOutputOptions & INCHI_OUT_STDINCHI))
    {
        bStdFormat = 0;
    }
    /* Some stereo */
    if (!(nMode & REQ_MODE_STEREO))
    {
        inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
        first = 0;
        inchi_ios_eprint(log_file, "  Stereo OFF\n");
    }
    else
    {
        if (!(TG_FLAG_POINTED_EDGE_STEREO & ip->bTautFlags))
        {
            if (first)
            {
                inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
                first = 0;
            }
            inchi_ios_eprint(log_file, "  Both ends of wedge point to stereocenters\n");
        }
    }
    if (ip->bDoNotAddH)
    {
        if (first)
        {
            inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
        }
        inchi_ios_eprint(log_file, "  Do not add H\n");
    }


#if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
    if (ip->bRenumber == 1)
    {
        inchi_ios_eprint(log_file, "\nGenerate InChI upon random atom renumbering\n\n");
    }
#endif

#if ( UNDERIVATIZE == 1 )
    if (ip->bUnderivatize == 1)
    {
        inchi_ios_eprint(log_file, "\nConvert input structure to derivative precursor before InChI calculation\n\n");
    }
    else if (ip->bUnderivatize == 3)
    {
        inchi_ios_eprint(log_file, "\nOutputs derivative information for the input structure\n\n");
    }
#endif



    /*  Generation/conversion indicator */
    if (bStdFormat)
    {
        if (!(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct)
        {
            inchi_ios_eprint(log_file, "Generating standard InChI\n");
        }

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) && !defined(TARGET_EXE_USING_API) )
        /* effective only in command line program InChI or stdInChI */
        else if (bInChI2Struct)
        {
            inchi_ios_eprint(log_file, "Converting InChI(s) to structure(s) in %s\n",
                (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) ?
                "MOL format" : "aux. info format");
        }
#endif
    }
    else
    {
        inchi_ios_eprint(log_file, "Generating non-standard InChI with the options: \n");
    }

    /* SDfile output */
    if (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)
    {
        inchi_ios_eprint(log_file,
            "Output SDfile only without stereochemical information and atom coordinates%s\n",
            (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT) ?
            "\n(write H isotopes as D, T)" : "");
    }

    /* Fixed/Mobile H */
    if (!bStdFormat)
    {
        if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_BASIC | REQ_MODE_TAUT))
        {
            inchi_ios_eprint(log_file, "  Mobile H Perception OFF (include FixedH layer)\n");
        }
        else if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_TAUT))
        {
            inchi_ios_eprint(log_file, "  Mobile H Perception ON  (omit FixedH layer)\n");
        }
        else if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_BASIC))
        {
            inchi_ios_eprint(log_file, "  Mobile H ignored\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "  Undefined Mobile H mode\n");
        }

        if ((ip->bTautFlags & TG_FLAG_VARIABLE_PROTONS))
        {
            if (!(ip->bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS))
            {
                inchi_ios_eprint(log_file, "  Disabled Aggressive (De)protonation\n");
            }
        }

#if ( FIND_RING_SYSTEMS != 1 )
        inchi_ios_eprint(log_file, "  %s5-, 6-, 7-memb. ring taut. ignored\n", i ? "; " : "");
#endif

        /* RecMet */
        if (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD)
        {
            if (ip->bTautFlags & TG_FLAG_RECONNECT_COORD)
            {
                inchi_ios_eprint(log_file, "  Include bonds to metals\n");
            }
            else
            {
                inchi_ios_eprint(log_file, "  Do not reconnect metals (omit RecMet layer)\n");
            }
        }
        else
        {
            inchi_ios_eprint(log_file, "  Do not disconnect metals\n");
        }

        /* Isotopic - always ON, output disabled. 09-17-2009*/
        /*
        if ( nMode & REQ_MODE_ISO )
        inchi_ios_eprint( log_file, "  Isotopic ON\n");
        else if ( nMode & REQ_MODE_NON_ISO )
        inchi_ios_eprint( log_file, "  Isotopic OFF\n");
        */

#if ( FIX_ADJ_RAD == 1 )
        if (ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS)
        {
            inchi_ios_eprint(log_file, "Fix Adjacent Radicals\n");
        }
#endif

        /* Stereo */
        if (nMode & REQ_MODE_STEREO)
        {
            inchi_ios_eprint(log_file, "  %s%s%s%sStereo ON\n",
                (nMode & REQ_MODE_NOEQ_STEREO) ? "Slow " : "",
                (nMode & REQ_MODE_REDNDNT_STEREO) ? "Redund. " : "",
                (nMode & REQ_MODE_NO_ALT_SBONDS) ? "No AltBond " : "",

                (nMode & REQ_MODE_RACEMIC_STEREO) ? "Racemic " :
                (nMode & REQ_MODE_RELATIVE_STEREO) ? "Relative " :
                (nMode & REQ_MODE_CHIR_FLG_STEREO) ? "Chiral Flag " : "Absolute ");

            if (0 == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
            {
                inchi_ios_eprint(log_file, "  Include undefined/unknown stereogenic centers and bonds\n");
            }
            else if (REQ_MODE_SC_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
            {
                inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic centers\n");
            }
            else if (REQ_MODE_SB_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
            {
                inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic bonds\n");
            }
            else
            {
                /*case REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
                inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic centers and bonds\n");
            }
            if (0 != (nMode & REQ_MODE_DIFF_UU_STEREO))
            {
                inchi_ios_eprint(log_file, "  Make labels for unknown and undefined stereo different\n");
            }

#if ( defined(MIN_SB_RING_SIZE) && MIN_SB_RING_SIZE > 0 )
            k = (ip->nMode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT;
            if (bRELEASE_VERSION != 1 || k != MIN_SB_RING_SIZE)
            {
                if (k >= 3)
                {
                    inchi_ios_eprint(log_file, "  Min. stereobond ring size: %d\n", k);
                }
                else
                {
                    inchi_ios_eprint(log_file, "  Min. stereobond ring size: NONE\n");
                }
            }
#endif
        }   /* Stereo */
    }   /* !bStdFormat */

    if (!bStdFormat)
    {
        if (TG_FLAG_KETO_ENOL_TAUT & ip->bTautFlags)
        {
            inchi_ios_eprint(log_file, "  Account for keto-enol tautomerism\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "  Do not account for keto-enol tautomerism\n");
        }
        if (TG_FLAG_1_5_TAUT & ip->bTautFlags)
        {
            inchi_ios_eprint(log_file, "  Account for 1,5-tautomerism\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "  Do not account for 1,5-tautomerism\n");
        }

        if (TG_FLAG_PT_22_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_22_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_22_00 tautomerism\n");
        if (TG_FLAG_PT_16_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_16_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_16_00 tautomerism\n");
        if (TG_FLAG_PT_06_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_06_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_06_00 tautomerism\n");
        if (TG_FLAG_PT_39_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_39_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_39_00 tautomerism\n");
        if (TG_FLAG_PT_13_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_13_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_13_00 tautomerism\n");
        if (TG_FLAG_PT_18_00 & ip->bTautFlags)
            inchi_ios_eprint(log_file, "  Account for PT_18_00 tautomerism\n");
        else
            inchi_ios_eprint(log_file, "  Do not account for PT_18_00 tautomerism\n");

        if (developer_options)
        {
            if (TG_FLAG_PHOSPHINE_STEREO & ip->bTautFlags)
            {
                inchi_ios_eprint(log_file, "  Include phosphine stereochemistry\n");
            }
            else
            {
                inchi_ios_eprint(log_file, "  Do not include phosphine stereochemistry\n");
            }
            if (TG_FLAG_ARSINE_STEREO & ip->bTautFlags)
            {
                inchi_ios_eprint(log_file, "  Include arsine stereochemistry\n");
            }
            else
            {
                inchi_ios_eprint(log_file, "  Do not include arsine stereochemistry\n");
            }
            if (!(TG_FLAG_FIX_SP3_BUG & ip->bTautFlags))
            {
                inchi_ios_eprint(log_file, "  Turned OFF fix of bug leading to missing or undefined sp3 parity\n");
            }
            if (!(TG_FLAG_FIX_ISO_FIXEDH_BUG & ip->bTautFlags))
            {
                inchi_ios_eprint(log_file, "  Turned OFF bug-fixes found after v.1.02b release\n");
            }
            if (!(ip->bFixNonUniformDraw))
            {
                inchi_ios_eprint(log_file, "  Turned OFF fixes of non-uniform drawing issues\n");
            }
            if (!(TG_FLAG_MOVE_POS_CHARGES & ip->bTautFlags))
            {
                inchi_ios_eprint(log_file, "  MovePos turned OFF\n");
            }
        }

    } /* !bStdFormat */

    if (ip->bCalcInChIHash != INCHIHASH_NONE)
    {
        if (bStdFormat)
        {
            inchi_ios_eprint(log_file, "Generating standard InChIKey\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "Generating InChIKey\n");
        }
        if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1)
        {
            inchi_ios_eprint(log_file, "Generating hash extension (1st block)\n");
        }
        else if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2)
        {
            inchi_ios_eprint(log_file, "Generating hash extension (2nd block)\n");
        }
        else if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2)
        {
            inchi_ios_eprint(log_file, "Generating hash extension (two blocks)\n");
        }
    }

    if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
    {
        inchi_ios_eprint(log_file, "Saving InChI creation options");
        if (bStdFormat)
        {
            inchi_ios_eprint(log_file, " suppressed for standard InChI");
            /* NB: actual suppression takes place on InChI serialization */
            /* (as on e.g. Inchi2Inchi conversion it may appear that we create non-std */
            /*  InChI instead of standard one) */
        }
        inchi_ios_eprint(log_file, "\n");
    }

    if (ip->bAllowEmptyStructure)
    {
        inchi_ios_eprint(log_file, "Issue warning on empty structure\n");
    }

    /* Input */
    if (ip->nInputType)
    {
        inchi_ios_eprint(log_file, "Input format: %s",
            ip->nInputType == INPUT_MOLFILE ? "MOLfile" :
            ip->nInputType == INPUT_SDFILE ? "SDfile" :
#if ( READ_INCHI_STRING == 1 )
            ip->nInputType == INPUT_INCHI ? "InChI (plain identifier)" :
#endif
            ip->nInputType == INPUT_INCHI_PLAIN ? "InChI AuxInfo (plain)" : "Unknown");
        if ((ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE) &&
            ip->bGetMolfileNumber)
        {
            inchi_ios_eprint(log_file, "  (attempting to read Molfile number)");
        }
        inchi_ios_eprint(log_file, "\n");
    }

    if (ip->szSdfDataHeader[0] && ip->nInputType != INPUT_SDFILE)
    {
        inchi_ios_eprint(log_file, "  SDfile data header: \"%s\"\n", ip->szSdfDataHeader);
    }

    /* Output */
    inchi_ios_eprint(log_file, "Output format: %s%s\n",
        (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) ? "Plain text" :

        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && bInChI2Struct) ? "SDfile only (without stereochemical info and atom coordinates)" :
        ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct) ? "SDfile only" : "Unknown",

        ((ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) &&
            (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)) ? ", tabbed" : "");

#if ( bRELEASE_VERSION == 1 )
    if (ip->bCtPredecessors || ip->bAbcNumbers)
    {
        if (ip->bCtPredecessors && ip->bAbcNumbers)
        {
            inchi_ios_eprint(log_file, "Representation: Compressed\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "Connection table: %s, %s\n",
                ip->bCtPredecessors ? "Predecessor_numbers(closures)" : "Canon_numbers(branching, ring closures)",
                ip->bAbcNumbers ? "Shorter alternative" : "Numerical");
        }
    }
#else
    if ((bRELEASE_VERSION != 1) || ip->bCtPredecessors || ip->bAbcNumbers)
    {
        inchi_ios_eprint(log_file, "Connection table: %s, %s\n",
            ip->bCtPredecessors ? "Predecessor_numbers(closures)" : "Canon_numbers(branching, ring closures)",
            ip->bAbcNumbers ? "Shorter alternative" : "Numerical");
    }
    else
    {
        inchi_ios_eprint(log_file, "Representation: Numerical");
    }
#endif

    if (ip->bNoWarnings)
    {
        inchi_ios_eprint(log_file, "Warnings suppressed\n");
    }

    if (ip->bHideInChI)
    {
        inchi_ios_eprint(log_file, "Printing InChI string itself suppressed\n");
    }
    if (ip->bMergeHash)
    {
        inchi_ios_eprint(log_file, "InChIKey combined with extra hash(es)\n");
    }



    if (!(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY))
    {
        if (ip->bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO)
        {
            inchi_ios_eprint(log_file, "Aux. info suppressed\n");
        }
        else if (ip->bINChIOutputOptions & INCHI_OUT_SHORT_AUX_INFO)
        {
            inchi_ios_eprint(log_file, "Minimal Aux. info\n");
        }
        else
        {
            inchi_ios_eprint(log_file, "Full Aux. info\n");
        }
    }
    if (ip->first_struct_number > 1)
    {
        inchi_ios_eprint(log_file, "Skipping %ld structure%s\n", ip->first_struct_number - 1, ip->first_struct_number == 2 ? "" : "s");
    }
    if (ip->last_struct_number > 0)
    {
        inchi_ios_eprint(log_file, "Terminate after structure #%ld\n", ip->last_struct_number);
    }
    if (ip->bSaveWarningStructsAsProblem && ip->path[3] && ip->path[3][0])
    {
        inchi_ios_eprint(log_file, "Saving warning structures into the problem file\n");
    }
    if (ip->bSaveAllGoodStructsAsProblem && ip->path[3] && ip->path[3][0])
    {
        inchi_ios_eprint(log_file, "Saving only all good structures into the problem file\n");
    }

    if (ip->bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR)
    {
        inchi_ios_eprint(log_file, "Print empty InChI if generation fails\n");
    }
    if (ip->bINChIOutputOptions2 & INCHI_OUT_MISMATCH_AS_ERROR)
    {
        inchi_ios_eprint(log_file, "Consider problem/mismatch on InChI conversion as error\n");
    }
    if (ip->msec_MaxTime)
    {
        inchi_ios_eprint(log_file, "Timeout per structure: %ld msec\n", ip->msec_MaxTime);
        /*
        unsigned long seconds = ip->msec_MaxTime / 1000;
        unsigned long milliseconds = (ip->msec_MaxTime%1000);
        inchi_ios_eprint( log_file, "Timeout per structure: %lu/*.%03lu sec\n", seconds, milliseconds); -- djb-rwth: ignoring LLVM warning
        inchi_ios_eprint( log_file, "Timeout per structure: %lu sec\n", seconds );
        */
    }
    else
    {
        inchi_ios_eprint(log_file, "No timeout\n");
    }

    if (ip->bLooseTSACheck)
    {
        inchi_ios_eprint(log_file, "Relax criteria of ambiguous drawing for in-ring stereo centers\n");
    }

    {
        int maxna = NORMALLY_ALLOWED_INP_MAX_ATOMS;
        if (ip->bLargeMolecules)
        {
            inchi_ios_eprint(log_file, "Experimental mode: ");
            maxna = MAX_ATOMS;
        }
        inchi_ios_eprint(log_file, "Up to %d atoms per structure\n", maxna);
    }

    if (ip->bPolymers != POLYMERS_NO)
    {
        inchi_ios_eprint(log_file, "Experimental mode: Treating polymers");

        if (ip->bPolymers == POLYMERS_MODERN)
        {
            ;
        }
        else if (ip->bPolymers == POLYMERS_LEGACY)
        {
            inchi_ios_eprint(log_file, " (v. 1.05 legacy mode)");
        }
        else if (ip->bPolymers == POLYMERS_LEGACY_PLUS)
        {
            inchi_ios_eprint(log_file, " (v. 1.05 legacy mode with senior link placed at start)");
        }
        if (ip->bFoldPolymerSRU)
        {
            inchi_ios_eprint(log_file, "; CRU folding enabled");
        }
    }
    inchi_ios_eprint(log_file, "\n");

    if (ip->bNPZz == 1)
    {
        inchi_ios_eprint(log_file, "Allowing non-polymer Zz pseudo atoms\n");
    }
    if (ip->bStereoAtZz == 1)
    {
        inchi_ios_eprint(log_file, "Allowing stereo at atoms connected to Zz\n");
    }

    /*  Report debug modes */
#if ( bRELEASE_VERSION != 1 )
    inchi_ios_eprint(log_file, "Release version = NO\n");
#endif
#if ( TRACE_MEMORY_LEAKS == 1 && defined(_DEBUG) )
    inchi_ios_eprint(log_file, "Tracing memory leaks (SLOW)\n");
#endif

#if (BUILD_WITH_ENG_OPTIONS==1)
    inchi_ios_eprint(log_file, "! Working in engineering mode\n");
#endif
    inchi_ios_eprint(log_file, "\n");

#if ( bRELEASE_VERSION != 1 )
#if ( FIND_RING_SYSTEMS == 1 )
    inchi_ios_eprint(log_file, "Find ring systems=Y\nTautomers:\n");
    inchi_ios_eprint(log_file, " 4-pyridinol=%s\n", TAUT_4PYRIDINOL_RINGS == 1 ? "Y" : "N");
    inchi_ios_eprint(log_file, " pyrazole=%s\n", TAUT_PYRAZOLE_RINGS == 1 ? "Y" : "N");
    inchi_ios_eprint(log_file, " tropolone=%s\n", TAUT_TROPOLONE_7 == 1 ? "Y" : "N");
    inchi_ios_eprint(log_file, " tropolone-5=%s\n", TAUT_TROPOLONE_5 == 1 ? "Y" : "N");
    inchi_ios_eprint(log_file, "Only chain attachments to tautomeric rings=%s\n", TAUT_RINGS_ATTACH_CHAIN == 1 ? "Y" : "N");
#endif

    if (ip->bGetSdfileId)
    {
        inchi_ios_eprint(log_file, "Extracting SDfile IDs\n");
    }
    inchi_ios_eprint(log_file, "\nDbg: MOVE_CHARGES=%d\n",
        0 != (ip->bTautFlags & TG_FLAG_MOVE_POS_CHARGES));
    inchi_ios_eprint(log_file, "     REPLACE_ALT_WITH_TAUT=%d; NEUTRALIZE_ENDPOINTS=%d; BNS_PROTECT_FROM_TAUT=%d\n",
        REPLACE_ALT_WITH_TAUT, NEUTRALIZE_ENDPOINTS, BNS_PROTECT_FROM_TAUT);
    inchi_ios_eprint(log_file, "     DISCONNECT_SALTS=%d;   TEST_TAUT_SALTS=%d;    TEST_TAUT2_SALTS=%d\n",
        0 != (ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS),
        0 != (ip->bTautFlags & TG_FLAG_TEST_TAUT__SALTS),
        0 != (ip->bTautFlags & TG_FLAG_TEST_TAUT2_SALTS));
    inchi_ios_eprint(log_file, "     CHARGED_ACID_TAUT_ONLY=%d MERGE_TAUT_SALTS=%d\n",
        0 == (ip->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O),
        0 != (ip->bTautFlags & TG_FLAG_MERGE_TAUT_SALTS));
    inchi_ios_eprint(log_file, "     DISCONNECT_COORD=%d\n",
        0 != (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD));
#endif /* ( bRELEASE_VERSION != 1 ) */

    return 0;
}


/****************************************************************************/
void HelpCommandLineParms(INCHI_IOSTREAM* f)
{
    if (!f)
    {
        return;
    }

#if ( bRELEASE_VERSION == 1 )
    inchi_ios_print_nodisplay(f,
#ifdef TARGET_EXE_USING_API
        "%s %-s\n%-s Build (%-s%-s) of %s %-s %-s\n\nUsage:\ninchi_main inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n",
        APP_DESCRIPTION, INCHI_SRC_REV,
        INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
        RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***",
        INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#else
        "%s %-s\n%-s Build (%-s%-s) of %s %-s %-s\n\nUsage:\ninchi-1 inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n",
        APP_DESCRIPTION, INCHI_SRC_REV,
        INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
        RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***",
        INCHI_OPTION_PREFX, INCHI_OPTION_PREFX );

#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay(f,
        "inchi-1 inputFiles... %cAMI [%coption[ %coption...]]\n",
        INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
#endif
#endif

    inchi_ios_print_nodisplay(f, "\nOptions:\n");

    inchi_ios_print_nodisplay(f, "\nInput\n");
    inchi_ios_print_nodisplay(f, "  STDIO       Use standard input/output streams\n");
    inchi_ios_print_nodisplay(f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    inchi_ios_print_nodisplay(f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");

    inchi_ios_print_nodisplay(f, "  START:n     Start at n-th input structure\n");
    inchi_ios_print_nodisplay(f, "  END:n       Stop after n-th input structure\n");
    inchi_ios_print_nodisplay(f, "  RECORD:n    Treat only n-th input structure\n");

#if ( BUILD_WITH_AMI == 1 )
    inchi_ios_print_nodisplay(f, "  AMI         Allow multiple input files (wildcards supported)\n");
    inchi_ios_print_nodisplay(f, "  AMIOutStd   Write output to stdout (in AMI mode)\n");
    inchi_ios_print_nodisplay(f, "  AMILogStd   Write log to stderr (in AMI mode)\n");
    inchi_ios_print_nodisplay(f, "  AMIPrbNone  Suppress creation of problem files (in AMI mode)\n");
#endif


    inchi_ios_print_nodisplay(f, "Output\n");
    inchi_ios_print_nodisplay(f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    inchi_ios_print_nodisplay(f, "  NoWarnings  Suppress all warning messages\n");
    inchi_ios_print_nodisplay(f, "  AuxNone     Omit auxiliary information\n");
    inchi_ios_print_nodisplay(f, "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n");
    inchi_ios_print_nodisplay(f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
    inchi_ios_print_nodisplay(f, "  MergeHash   Combine InChIKey with extra hash(es) if present\n");
    inchi_ios_print_nodisplay(f, "  NoInChI     Do not print InChI string itself\n");
#ifndef TARGET_EXE_USING_API
    inchi_ios_print_nodisplay(f, "  OutErrInChI On fail, print empty InChI (default: nothing)\n");
#endif
#if ( defined(_WIN32) && !defined(COMPILE_ANSI_ONLY) && !defined(TARGET_API_LIB) ) /* djb-rwth: check if this is working on GCC for Windows */
    inchi_ios_print_nodisplay(f, "  D           Display the structures\n");
    inchi_ios_print_nodisplay(f, "  EQU         Display sets of identical components\n");
    inchi_ios_print_nodisplay(f, "  Fnumber     Set display Font size in number of points\n");
#endif
    inchi_ios_print_nodisplay(f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
#if ( SDF_OUTPUT_DT == 1 )
    inchi_ios_print_nodisplay(f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
#endif

    inchi_ios_print_nodisplay(f, "Structure perception\n");
    inchi_ios_print_nodisplay(f, "  SNon        Exclude stereo (default: include absolute stereo)\n");
    inchi_ios_print_nodisplay(f, "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n");
    inchi_ios_print_nodisplay(f, "  LooseTSACheck   Relax criteria of ambiguous drawing for in-ring tetrahedral stereo\n");
    inchi_ios_print_nodisplay(f, "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n");
#ifndef USE_STDINCHI_API
    inchi_ios_print_nodisplay(f, "Stereo perception modifiers (non-standard InChI)\n");
    inchi_ios_print_nodisplay(f, "  SRel        Relative stereo\n");
    inchi_ios_print_nodisplay(f, "  SRac        Racemic stereo\n");
    inchi_ios_print_nodisplay(f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n");

    inchi_ios_print_nodisplay(f, "Customizing InChI creation (non-standard InChI)\n");
    inchi_ios_print_nodisplay(f, "  SUU         Always include omitted unknown/undefined stereo\n");
    inchi_ios_print_nodisplay(f, "  SLUUD       Make labels for unknown and undefined stereo different\n");
    inchi_ios_print_nodisplay(f, "  RecMet      Include reconnected metals results\n");
    inchi_ios_print_nodisplay(f, "  FixedH      Include Fixed H layer\n");

    inchi_ios_print_nodisplay(f, "  KET         Consider keto-enol tautomerism (experimental)\n");
    inchi_ios_print_nodisplay(f, "  15T         Consider 1,5-tautomerism (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_06_00    Consider 1,3 heteroatom shift (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_13_00    Consider keten-ynol exchange (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_16_00    Consider nitroso-oxime tautomerism (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_18_00    Consider cyanic/iso-cyanic acids (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_22_00    Consider imine/imine tautomerism (experimental)\n");
    inchi_ios_print_nodisplay(f, "  PT_39_00    Consider nitrone/azoxy or Behrend rearrangement (experimental)\n");


    inchi_ios_print_nodisplay(f, "Generation\n");
    inchi_ios_print_nodisplay(f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    inchi_ios_print_nodisplay(f, "  WMnumber    Set time-out per structure in milliseconds (int); WM0 means unlimited\n");
    inchi_ios_print_nodisplay(f, "  LargeMolecules Treat molecules up to 32766 atoms (experimental)\n");
    inchi_ios_print_nodisplay(f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);
    /*inchi_ios_print_nodisplay( f, "  MismatchIsError Treat problem/mismatch on inchi2struct conversion as error\n");*/

    inchi_ios_print_nodisplay(f, "  Polymers    Allow processing of polymers (experimental)\n");
    inchi_ios_print_nodisplay(f, "  Polymers105 Allow processing of polymers (experimental, legacy mode of v. 1.05)\n");
    inchi_ios_print_nodisplay(f, "  FoldCRU     Fold polymer CRU if inner repeats occur\n");
    inchi_ios_print_nodisplay(f, "  NoFrameShift Disable polymer CRU frame shift\n");
    inchi_ios_print_nodisplay(f, "  NoEdits     Disable polymer CRU frame shift and folding\n");
    inchi_ios_print_nodisplay(f, "  NPZz        Allow non-polymer-related Zz atoms (pseudo element placeholders)\n");
    inchi_ios_print_nodisplay(f, "  SAtZz       Allow stereo at atoms connected to Zz(default: disabled)\n");
#endif


    inchi_ios_print_nodisplay(f, "  Key         Generate InChIKey\n");
    inchi_ios_print_nodisplay(f, "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n");
    inchi_ios_print_nodisplay(f, "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n");

    inchi_ios_print_nodisplay(f, "Conversion\n");
#ifdef TARGET_EXE_USING_API
    inchi_ios_print_nodisplay(f, "  InChI2Struct Test mode: Mol/SDfile -> %s -> Structure -> (%s+AuxInfo)\n", INCHI_NAME, INCHI_NAME);
#else
    inchi_ios_print_nodisplay(f, "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n");
    inchi_ios_print_nodisplay(f, "  InChI2InChI  Convert  Convert %s string(s) into %s string(s)\n", INCHI_NAME, INCHI_NAME);
#endif

#if (BUILD_WITH_ENG_OPTIONS==1)
    inchi_ios_print_nodisplay(f, "Engineering/hidden\n");
#ifdef TARGET_EXE_USING_API
    inchi_ios_print_nodisplay(f, "  InChI2InChI  Test mode: Mol/SDfile -> %s -> %s\n", INCHI_NAME, INCHI_NAME);
#endif

    /*inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n"); */
    /*inchi_ios_print_nodisplay( f, "    FULL        Standard set of options for Full Verbose Output\n");*/
    /*inchi_ios_print_nodisplay( f, "    MIN         Standard set of options for Minimal Concise Output\n");*/

#if ALLOW_SUBSTRUCTURE_FILTERING==1
    inchi_ios_print_nodisplay(f, "  FilterSS    Select input SDF records using (hard-coded) substructure filter\n");
    inchi_ios_print_nodisplay(f, "  InvFilterSS Invert match for (hard-coded) substructure filter\n");
#endif

    inchi_ios_print_nodisplay(f, "  Compress    Compressed output\n");
    inchi_ios_print_nodisplay(f, "  MERGE       Use bMergeAllInputStructures\n");
    inchi_ios_print_nodisplay(f, "  PGO         Use bSaveAllGoodStructsAsProblem\n");
    inchi_ios_print_nodisplay(f, "  DCR         Use bDisplayCompositeResults\n");
    inchi_ios_print_nodisplay(f, "  DSB         Use REQ_MODE_NO_ALT_SBONDS \n");
    inchi_ios_print_nodisplay(f, "  NOHDR       Use bNoStructLabels\n");
    inchi_ios_print_nodisplay(f, "  NoVarH      Set bTgFlagVariableProtons=0\n");
    inchi_ios_print_nodisplay(f, "  NOUUSB      Use REQ_MODE_SB_IGN_ALL_UU\n");
    inchi_ios_print_nodisplay(f, "  NOUUSC      Use REQ_MODE_SC_IGN_ALL_UU\n");
    inchi_ios_print_nodisplay(f, "  FixRad      Set bFixAdjacentRad\n");
    inchi_ios_print_nodisplay(f, "  TestRenum   Generate InChI upon random atom renumbering\n");
    inchi_ios_print_nodisplay(f, "  DoDRV       Set bUnderivatize=1\n");
    inchi_ios_print_nodisplay(f, "  DoDrvReport Set bUnderivatize=3\n");
    inchi_ios_print_nodisplay(f, "  DoR2C       Set bRing2Chain\n");
    inchi_ios_print_nodisplay(f, "  DoneOnly    Set bIgnoreUnchanged\n");
    inchi_ios_print_nodisplay(f, "  NoADP       Set bTgFlagHardAddRenProtons=0\n");
    inchi_ios_print_nodisplay(f, "  MOVEPOS:0|1 Set bMovePositiveCharges\n");
    inchi_ios_print_nodisplay(f, "  RSB:n       Set nMinDbRinSize\n");
    inchi_ios_print_nodisplay(f, "  DISCONSALT:0|1     Set bDisconnectSalts\n");
    inchi_ios_print_nodisplay(f, "  DISCONMETAL:0|1    Set bDisconnectCoord\n");
    inchi_ios_print_nodisplay(f, "  DISCONMETALCHKVAL:0|1 Set bDisconnectCoordChkVal \n");
    inchi_ios_print_nodisplay(f, "  RECONMETAL:0|1     Set bReconnectCoord\n");
    inchi_ios_print_nodisplay(f, "  MERGESALTTG:0|1    Set bMergeSaltTGroups\n");
    inchi_ios_print_nodisplay(f, "  UNCHARGEDACIDS:0|1 Set bUnchargedAcidTaut \n");
    inchi_ios_print_nodisplay(f, "  ACIDTAUT:0|1|2     Set bAcidTautomerism\n");
    inchi_ios_print_nodisplay(f, "  AUXINFO:0|1|2      Set AuxInfo print options\n");
    inchi_ios_print_nodisplay(f, "  KeepBalanceP...  \n");
    inchi_ios_print_nodisplay(f, "  SDFID       ...\n");
    inchi_ios_print_nodisplay(f, "  PLAINP      ....\n");
    inchi_ios_print_nodisplay(f, "  ANNPLAIN    ....\n");

#endif

#endif
}

#define fprintf2 inchi_fprintf

#ifndef TARGET_API_LIB
/****************************************************************************/
int OpenFiles(FILE** inp_file,
    FILE** out_file,
    FILE** log_file,
    FILE** prb_file,
    INPUT_PARMS* ip)
{
    /*
    -- Files --
    ip->path[0] => Input
    ip->path[1] => Output (INChI)
    ip->path[2] => Log
    ip->path[3] => Problem structures
    ip->path[4] => Errors file (ACD Labs)

    */

    /*  Logfile (open as early as possible) */
    if (!ip->path[2] || !ip->path[2][0])
    {
        fprintf2(stderr, "%s %-s\n%-s Build (%-s%-s) of %s %s %-s\n\n",
            APP_DESCRIPTION, INCHI_SRC_REV,
            INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
            RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***");
        fprintf2(stderr, "Log file not specified. Using standard error output.\n");
        *log_file = stderr;
    }
    else if (!(*log_file = fopen(ip->path[2], "w")))
    {
        fprintf2(stderr, "%s %-s\n%-s Build (%-s%-s) of %s %s%-s\n\n",
            APP_DESCRIPTION, INCHI_SRC_REV,
            INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
            RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***");
        fprintf2(stderr, "Cannot open log file '%s'. Using standard error output.\n", ip->path[2]);
        *log_file = stderr;
    }
    else
    {
        fprintf2(*log_file, "%s %-s\n%-s Build (%-s%-s) of %s %s%-s\n\n",
            APP_DESCRIPTION, INCHI_SRC_REV,
            INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
            RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***");
        fprintf2(*log_file, "Opened log file '%s'\n", ip->path[2]);
    }

    /* Input file */
    if ((ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE ||
        ip->nInputType == INPUT_INCHI || ip->nInputType == INPUT_INCHI_PLAIN)
        &&
        ip->num_paths > 0)
    {
        const char* fmode = NULL;

#if ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        /* compilers that definitely allow fopen "rb" (binary read) mode */
        fmode = "rb";
        if (!ip->path[0] || !ip->path[0][0] || !(*inp_file = fopen(ip->path[0], "rb")))
        {
            fprintf2(*log_file, "Cannot open input file '%s'. Terminating.\n", ip->path[0] ? ip->path[0] : "<No name>");
            goto exit_function;
        }
        else
        {
            fprintf2(*log_file, "Opened input file '%s'\n", ip->path[0]);
        }

#else

        if (!ip->path[0] || !ip->path[0][0] || !(*inp_file = fopen(ip->path[0], "r")))
        {
            fprintf2(*log_file, "Cannot open input file '%s'. Terminating.\n", ip->path[0] ? ip->path[0] : "<No Name>");
            goto exit_function;
        }
        else
        {
            fprintf2(*log_file, "Opened input file '%s'\n", ip->path[0]);
        }
        fmode = "r";
#endif /* ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) ) */

        DetectInputINChIFileType(inp_file, ip, fmode);
    }

    else if ((ip->nInputType != INPUT_MOLFILE &&
        ip->nInputType != INPUT_SDFILE &&
        ip->nInputType != INPUT_INCHI &&
        /* post-1.02b */
        ip->nInputType != INPUT_INCHI_PLAIN))
    {
        fprintf2(*log_file, "Input file type not specified. Terminating.\n");
        goto exit_function;
    }
    else
    {
        fprintf2(*log_file, "Input file not specified. Using standard input.\n");
        *inp_file = stdin;
    }

    /*  Output file */
    if (!ip->path[1] || !ip->path[1][0])
    {
        fprintf2(*log_file, "Output file not specified. Using standard output.\n");
        *out_file = stdout;
    }
    else
    {
        if (!(*out_file = fopen(ip->path[1], "w")))
        {
            fprintf2(*log_file, "Cannot open output file '%s'. Terminating.\n", ip->path[1]);
            goto exit_function;
        }
        else
        {
            fprintf2(*log_file, "Opened output file '%s'\n", ip->path[1]);
            if ((ip->bINChIOutputOptions & (INCHI_OUT_PLAIN_TEXT)) &&
                *inp_file != stdin &&
                !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) &&
                !ip->bNoStructLabels &&
                !(ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT))
            {
                PrintFileName("* Input_File: \"%s\"\n", *out_file, ip->path[0]);
            }
        }
    }

    /*  Problem file */
    if (ip->path[3] && ip->path[3][0])
    {
        const char* fmode = "w";

#if ( defined(_MSC_VER)&&defined(_WIN32) || defined(__BORLANDC__)&&defined(__WIN32__) || defined(__GNUC__)&&defined(__MINGW32__)&&defined(_WIN32) )
        fmode = "wb";
#endif

        if (!(*prb_file = fopen(ip->path[3], fmode)))
        {
            fprintf2(*log_file, "Cannot open problem file '%s'. Terminating.\n", ip->path[3]);
            goto exit_function;
        }
        else
        {
            fprintf2(*log_file, "Opened problem file '%s'\n", ip->path[3]);
        }
    }

    /*  Success */
    return 1;

exit_function:

    /*  Failed */
    return 0;
}


#define NUM_VERSIONS 7
#define LEN_VERSIONS 64


static int bMatchOnePrefix(int len, char* str, int lenPrefix[],
    char strPrefix[][LEN_VERSIONS], int numPrefix);


/****************************************************************************/
static int bMatchOnePrefix(int len, char* str,
    int lenPrefix[],
    char strPrefix[][LEN_VERSIONS],
    int numPrefix)
{
    int i;
    for (i = 0; i < numPrefix; i++)
    {
        if (len >= lenPrefix[i] &&
            !memcmp(str, strPrefix[i], lenPrefix[i]))
        {
            return 1;
        }
    }

    return 0;
}

/****************************************************************************/
int DetectInputINChIFileType(FILE** inp_file,
    INPUT_PARMS* ip,
    const char* fmode)
{
    char szLine[256], ret = 0;
    static char szPlnVersion[NUM_VERSIONS][LEN_VERSIONS]; /* = "INChI:1.1Beta/";*/
    static int  lenPlnVersion[NUM_VERSIONS];
    static char szPlnAuxVer[NUM_VERSIONS][LEN_VERSIONS]; /* = "AuxInfo:1.1Beta/";*/
    static int  lenPlnAuxVer[NUM_VERSIONS];
    static int  bInitialized = 0;
    int  bINChI_plain = 0, len, i;


    if (ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_INCHI)
    {
        return 1;
    }

    if (!bInitialized)
    {
        lenPlnVersion[0] = sprintf(szPlnVersion[0], "%s=%s/", INCHI_NAME, INCHI_VERSION);
        lenPlnVersion[1] = sprintf(szPlnVersion[1], "INChI=1.12Beta/");
        lenPlnVersion[2] = sprintf(szPlnVersion[2], "INChI=1.0RC/");
        lenPlnVersion[3] = sprintf(szPlnVersion[3], "InChI=1.0RC/");
        lenPlnVersion[4] = sprintf(szPlnVersion[4], "InChI=1/");
        lenPlnVersion[5] = sprintf(szPlnVersion[5], "MoChI=1a/");
        lenPlnVersion[6] = sprintf(szPlnVersion[6], "InChI=1S/");
        lenPlnAuxVer[0] = sprintf(szPlnAuxVer[0], "AuxInfo=%s/", INCHI_VERSION);
        lenPlnAuxVer[1] = sprintf(szPlnAuxVer[1], "AuxInfo=1.12Beta/");
        lenPlnAuxVer[2] = sprintf(szPlnAuxVer[2], "AuxInfo=1.0RC/");
        lenPlnAuxVer[3] = sprintf(szPlnAuxVer[3], "AuxInfo=1.0RC/");
        lenPlnAuxVer[4] = sprintf(szPlnAuxVer[4], "AuxInfo=1/");
        lenPlnAuxVer[5] = sprintf(szPlnAuxVer[5], "AuxInfo=1a/");
        lenPlnAuxVer[6] = sprintf(szPlnAuxVer[6], "AuxInfo=1/");
#if ( FIX_DALKE_BUGS == 1 )
        bInitialized = 1;
#endif
    }

    for (i = 0; i < 4; i++)
    {
        len = inchi_fgetsLfTab(szLine, sizeof(szLine) - 1, *inp_file);
        if (len < 0)
        {
            break;
        }
        if (bMatchOnePrefix(len, szLine, lenPlnVersion, szPlnVersion, NUM_VERSIONS) ||
            bMatchOnePrefix(len, szLine, lenPlnAuxVer, szPlnAuxVer, NUM_VERSIONS))
        {
            bINChI_plain++;
        }
    }

    if (bINChI_plain >= 2)
    {
        ip->nInputType = INPUT_INCHI_PLAIN;
        ret = 1;
    }

    fclose(*inp_file);
    *inp_file = fopen(ip->path[0], fmode);

    return ret;
}
#undef NUM_VERSIONS
#undef LEN_VERSIONS

#endif /* TARGET_API_LIB */
