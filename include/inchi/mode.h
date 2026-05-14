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


#ifndef _MODE_H_
#define _MODE_H_

#include <stdio.h>
#include <stdlib.h>
#include "bcf_s.h"


/*******************/
/*                 */
/*  BUILD TARGETS  */
/*                 */
/*******************/

/*
    One of targets below should be explicitly
    indicated to compiler (makefile/directive):

        TARGET_EXE_STANDALONE    Stand-alone executable inchi-1[.exe]
        TARGET_API_LIB           Library (libinchi) for using InChI API
                                 described in inchi_api.h
        TARGET_EXE_USING_API     Executable (INCHI_MAIN) which uses
                                 API library (e.g., libinchi.dll)
        TARGET_LIB_FOR_WINCHI    library for wInChI
*/


#if defined(TARGET_EXE_STANDALONE)
/* Uncomment the next line to perform atom ordering/renumbering tests (internal only) */
#define RENUMBER_ATOMS_AND_RECALC_V106 1
#ifdef      RENUMBER_ATOMS_AND_RECALC_V106
/* Comment the next line to print all renumbering changing the initial InChIKey */
#define STOP_AFTER_FIRST_CHANGE_ON_RENUMBERING 1
/* djb-rwth: adding full version number in the output -- GH issue #61 */
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (inchi-1 executable)"
/*#define APP_DESCRIPTION "InChI version 1, Software v. 1.06-PT6 (inchi-1 executable) \n*** INTERNAL TEST MODE: ATOM RENUMBERING TEST IS ACTIVE ***"*/
#else
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (inchi-1 executable)***"
/*#define APP_DESCRIPTION "InChI version 1, Software v. 1.06-PT6 (inchi-1 executable) \n*** UNOFFICIAL TEST VERSION: 6 PT TAUTO RULES AVAILABLE ***"*/
#endif

#elif defined(TARGET_API_LIB)
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (API Library)"

#elif defined(TARGET_EXE_USING_API)
#ifndef APP_DESCRIPTION
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (executable calling API Library)"
#endif

#elif defined(TARGET_LIB_FOR_WINCHI)
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (Library for wInChI GUI executable)"

#elif defined(TARGET_WINCHI)
#define APP_DESCRIPTION "InChI version 1, Software " CURRENT_VER " (wInChI GUI executable)"

#else
#error  No build target #defined, pls check compiler options... (TARGET_EXE_STANDALONE|TARGET_API_LIB|TARGET_EXE_USING_API|TARGET_LIB_FOR_WINCHI)
#endif


#ifndef TARGET_PLATFORM
#if defined(_WIN32)
#define TARGET_PLATFORM "Windows"
#else
#define TARGET_PLATFORM "Linux"
#endif
#endif


/****************************/
/*                          */
/*  BUILD OPTIONS/FEATURES  */
/*                          */
/****************************/

/*    Possible options are:

BUILD_LINK_AS_DLL
    Link library as a Win32 DLL or to eliminate inchi_stricmp duplication
    (use with TARGET_API_LIB or TARGET_EXE_USING_API)

BUILD_WITH_ENG_OPTIONS
    Expose engineering options

BUILD_WITH_AMI
    Turns on AMI (Allow Multiple Inputs) mode for standalone executable

    Select and uncomment whichever are necessary from the list below. */


/* #define BUILD_LINK_AS_DLL */


#ifndef BUILD_WITH_ENG_OPTIONS
#define BUILD_WITH_ENG_OPTIONS 0
#endif



#ifndef BUILD_WITH_AMI
/* this allows BUILD_WITH_AMI be #defined in a makefile */
#define BUILD_WITH_AMI 1
#endif
/* NB:  AMI mode is only for stand-alone executable */
#ifndef TARGET_EXE_STANDALONE
#ifdef BUILD_WITH_AMI
#undef BUILD_WITH_AMI
#endif
#endif

/* Smarter AMI for Windows */
/* Thanks, DT (2013-12-18) */
#if( BUILD_WITH_AMI == 1 )
#if( defined( _MSC_VER ) )
#define MSC_AMI 1
/* use MSC _findfirst(...), etc. Do not link with setargv.obj */
#endif
#if (BUILD_WITH_ENG_OPTIONS==1)
#define OUTPUT_FILE_EXT 1           /* options "/.": /.extOut /.extLog /.extPrb replace input file extension instead of adding .txt, .log, .prb 2013-12-18 DCh */
#define ALLOW_EMPTY_PATHS 1         /* let command line argument "" be interpreted as path */

#define DEL_EMPTY_OUTPUT  1         /* delete output file, problem file if the file is empty -- needed for /DoDrv /DoneOnly or to get rid of empty prb files */
#define SDF_OUTPUT_HETERO_VALENCE 1 /* output all hetero atom valences to SDF -- NIST output specific */

#endif
#endif

/* CML input is not supported started from v. 1.04 */
/* set ADD_CMLPPP to zero to override possble makefile define */
#define ADD_CMLPP 0


/*****************************/
/*                           */
/* COMPILE OPTIONS/FEATURES  */
/*                           */
/*****************************/

/*    Possible options are:

COMPILE_ANSI_ONLY
    Unconditionally force ANSI-89 C, no Win32 specific code

COMPILE_ADD_NON_ANSI_FUNCTIONS
    Use with COMPILE_ANSI_ONLY to add inchi_stricmp(), etc., see util.c

COMPILE_ALL_CPP
    allow C++ compilation/linkage of functions prototyped in .h files

MS VC compiler pragmas

    Select and uncomment whichever are necessary from the list below. */


/* #define COMPILE_ANSI_ONLY */
#if ( !defined(_MSC_VER) || defined(TARGET_API_LIB) || defined(TARGET_EXE_USING_API))  /* non-Microsoft GNU C, BCC, etc. compilers */
#ifndef COMPILE_ANSI_ONLY
#define COMPILE_ANSI_ONLY
#endif
#endif
#ifdef COMPILE_ANSI_ONLY
/*#define COMPILE_ADD_NON_ANSI_FUNCTIONS */
#endif

/*
#define COMPILE_ALL_CPP 1
*/

#ifdef _MSC_VER
/*
========== disable MS VC++ 6.0 Level 4 compiler warnings: ==============
 C4706: assignment within conditional expression
 C4127: conditional expression is constant
 C4244: '=' : conversion from 'int ' to '???', possible loss of data
 C4267: '=' : conversion from 'size_t' to 'int', possible loss of data
 C4701: local variable '???' may be used without having been initialized (removed)
 C4514: unreferenced inline/local function has been removed (C++)
 C4100: 'identifier' : unreferenced formal parameter
 C4786: 'identifier' : identifier was truncated to 'number' characters in the debug information
 C4996: 'identifier' was declared deprecated
========================================================================
*/
#pragma warning( disable : 4706 4127 4514 4100 4786 4996 4244 4267 )
#endif


/* To use allocator, uncomment the next line or use /D "USE_TBB_MALLOC=1" compile directive */
/*#define USE_TBB_MALLOC 1*/
#ifdef USE_TBB_MALLOC
#if !defined(TARGET_EXE_STANDALONE) && !defined(TARGET_LIB_FOR_WINCHI) && defined(_WIN32)
/* Use dynamic memory allocator from  Intel(R) Threading Building Blocks                    */
/* https://www.threadingbuildingblocks.org/                                                 */
/* see the header, tbbmalloc_proxy.h                                                        */
/* NB: under MS Visual Studio, add paths to platform-dependent TBB 32/64 bit libs          */
/* (tbbmalloc_proxy_debug.lib tbbmalloc_proxy.lib)                                          */
/* also ensure that final binaries access platform-dependent TBB 32/64 bit dlls             */
/* (tbbmalloc.dll tbbmalloc_proxy.dll) and MSVCR120.dll or other corresp. to TBB compile    */
/*                                                                                          */
#include "../../INCHI_API/tbb/tbbmalloc_proxy_for_inchi.h"
#endif
#endif


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/*********************/
/*                   */
/*  INCHI ALGORITHM  */
/*                   */
/*********************/

#define INCHI_VERSION       "1"
#define INCHI_NAME          "InChI"
#define INCHI_NAM_VER_DELIM "="

#ifdef _WIN32
#define   INCHI_OPTION_PREFX  '/'
#define   INCHI_PATH_DELIM    '\\'
#else
#define   INCHI_OPTION_PREFX  '-'
#define   INCHI_PATH_DELIM    '/'
#endif

#define   INCHI_ALT_OPT_PREFIX  '-'
#define   INCHI_ACD_LABS_PREFIX '-'

#define bRELEASE_VERSION  1    /* 1=> release version; comment out to disable */
#ifndef bRELEASE_VERSION
#define bRELEASE_VERSION  0    /* 0=> debug version */
#endif

/*#define RELEASE_IS_FINAL  0*/ /* 1=> pre-release version; comment out to disable */
#ifndef RELEASE_IS_FINAL
#define RELEASE_IS_FINAL  1    /* final release */
#endif

/* display (non-canonical) c-groups, display orig at numbers */
#if ( bRELEASE_VERSION == 1 )
#define DISPLAY_DEBUG_DATA_C_POINT 0  /* disabled release version for now */
#define DISPLAY_ORIG_AT_NUMBERS    1  /* 1 => in an uncanonicalized components display orig. atom numbers (default) */
#else
#define DISPLAY_DEBUG_DATA_C_POINT 1  /* debug: 1=>display (non-canonically numbered) c-groups, 0=>do not display */
#define DISPLAY_ORIG_AT_NUMBERS    1  /* 0 => in an uncanonicalized components display ordering atom numbers (debug) */
#endif

#if ( DISPLAY_DEBUG_DATA_C_POINT > 0 )
#define DISPLAY_DEBUG_DATA         DISPLAY_DEBUG_DATA_C_POINT
#endif


#define DISPLAY_ZZ_AS_STAR 1


/* BUG FIXES */

/**************************/
/* bug fixes in v1.00     */
/**************************/
#define FIX_ChCh_STEREO_CANON_BUG     1  /* 1=> (NEEDED) */
#define ADD_ChCh_STEREO_CANON_CHK     0  /* 1 is NOT needed; let it always be 0 */
#define FIX_ChCh_CONSTIT_CANON_BUG    1  /* 1=> (NEEDED) */
#define FIX_EITHER_STEREO_IN_AUX_INFO 1  /* 1=> fix bug: Either stereobond direction in Aux_Info; 0=> do not fix */
#define FIX_NORM_BUG_ADD_ION_PAIR     1  /* 1=> (NEEDED) fix bug: Miscount number of charges when creating an ion pair */
#define FIX_REM_PROTON_COUNT_BUG      1  /* 1=> (NEEDED) check for number of actually removed protons and issue an error if mismatch */
#define FIX_READ_AUX_MEM_LEAK         1
#define FIX_READ_LONG_LINE_BUG        1  /* 1=> (NEEDED) prevent failure when reading AuxInfo and InChI is too long */
#define FIX_N_V_METAL_BONDS_GPF       1  /* 1=> (NEEDED) InChI v1 GPF bug fix */
#define BNS_RAD_SEARCH                1  /* 1=> prevent normalization failures due to radical centers */

/*******************************/
/* bug fixes in post-v1.00     */
/*******************************/
#define FIX_ODD_THINGS_REM_Plus_BUG   0
#define FIX_N_MINUS_NORN_BUG          0
#define FIX_CANCEL_CHARGE_COUNT_BUG   0
#define FIX_2D_STEREO_BORDER_CASE     0
#define FIX_REM_ION_PAIRS_Si_BUG      0
#define FIX_STEREO_SCALING_BUG        0
#define FIX_EMPTY_LAYER_BUG           0
#define FIX_EITHER_DB_AS_NONSTEREO    0
#define FIX_BOND23_IN_TAUT            0
#define FIX_TACN_POSSIBLE_BUG         0
#define FIX_KEEP_H_ON_NH_ANION        0
#define FIX_AVOID_ADP                 0
/* may change InChI */
#define FIX_NUM_TG                    0  /* increase number of t-groups for isothiocyanate */
/* changes InChI for isothiocyanate */
#define FIX_CPOINT_BOND_CAP2          0

/*******************************/
/* bug fixes in post-v1.02b    */
/*******************************/

#define FIX_ISO_FIXEDH_BUG            1 /* (2007-09-24) 1=> Fix bug: missing fixed-H iso segment in case of single removed D(+) */
#define FIX_ISO_FIXEDH_BUG_READ       0 /* (2007-09-24) 1=> Accommodate this InChI bug in reading InChI */
#define FIX_DALKE_BUGS                1
#define FIX_TRANSPOSITION_CHARGE_BUG  1 /* (2008-01-02) fix bug that leads to missed charge in some cases when /o is present */
#define FIX_I2I_STEREOCONVERSION_BUG  1 /* (2008-03-06)   1=> Fix bug of i2i conversion SAbs-->(SRel||Srac) */
#define FIX_I2I_STEREOCONVERSION_BUG2 1 /* (2008-04-02)   1=> Fix bug of i2i conversion (missed empty /t) */
#define FIX_I2I_STEREOCONVERSION_BUG3 1 /* (2008-04-10)   1=> Fix bug of i2i conversion */
                                        /* (missed repeating /s in FI after F for multi-component case) */
#define FIX_TERM_H_CHRG_BUG           1 /* (2008-06-06) IPl) */
                                        /* fix bug: in some cases (dependent on ordering
                                        numbers), moving a charge from terminal H to heavy
                                        atom resulted in neutralizing H but not adjusting
                                        charge of heavy atom */

#define FIX_AROM_RADICAL              1 /* (2011-05-09) 1=> Fix bug which leads for different InChI */
                                        /* on atomic permitations for systems containing radical at */
                                        /* atom in aromatic ring */


/* Software version 1.07 */

/* INTENTIONALLY DISABLE 1.06 FIXES */
/*#define DISABLE_106_FIXES 1  */

#ifndef DISABLE_106_FIXES

#define FIX_STEREOCOUNT_ERR           1 /* (2018-01-09) Supplied by DT                              */
                                        /* Fix for InChI Error -30010 (STEREOCOUNT_ERR)             */
                                        /* appeared on PubChem CIDs 124897603, 124921144            */
                                        /* "The failure occurs when one of two or more              */
                                        /* constitutionally equivalent undefined stereocenters has  */
                                        /* been removed due to stereo equivalence of its two        */
                                        /* attachments."                                            */

#define FIX_IMPOSSIBLE_H_ISOTOPE_BUG  1 /* (2018-01-25) Bug reported by Andrew Dalke, inchi-discuss */
                                        /* Unrealistic H with mass difference 30 was consumed by    */
                                        /* InChI and sometimes resulted in memory corruption on     */
                                        /* accessing array of len == NUM_H_ISOTOPES that is 3.      */

#define FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG 1 
                                        /* (2018-05-03) thanks, DT; fix for 'olean spiro chirality' */


#define FIX_RENUM_BUG_FOR_CASE_OF_ACIDIC_OH_AT_P_PLUS 1
                                        /* Fix renumbering instability for case of acidic group     */
                                        /* near charged P and some other heteroatoms                */





/* Fixing issues reported by "Google Autofuzz project" and CURE53                                   */

/* internal v. 1.051, September 2019                                                                */
#define CHECK_STRTOL_ATNUMB 1
#define DISABLE_READ_COMPRESSED_INCHI 1
#define FIX_GAF_2019_1 1
#define FIX_GAF_2019_2 1

/* v. 1.06 October 2023 */
/*  Sep-Dec 2020 Fixed oss-fuzz issues:
    25604 25607 25609 25618 25726 25727 25728 25730 25731 25735
    25736 25737 25741 25830 25835 26514 27901 27903 27905
    Thanks to google/oss-fuzz for detecting and reporting the issues
 */

#define FIX_GAF_2020_GENERIC 1
#define FIX_OSS_FUZZ_25604 1
#define FIX_GAF_2020_25607 1
#define FIX_GAF_2020_25726 1
#define FIX_GAF_2020_25741 1

/* v. 1.061 August 2021 */
#define FIX_OSS_FUZZ_30162_30343 1
#define FIX_OSS_FUZZ_25734_28139 1

/* internal v. 1.052, January 2020 */
/* Thanks to CURE53 for detecting and reporting the issues */
#define FIX_CURE53_ISSUE_OOB_ALREADY_HAVE_THIS_MESSAGE 1
#define FIX_CURE53_ISSUE_HEAP_BUFFER_OVERFLOW_INCHITOINPATOM 1
#define FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO 1
                                        /* NB: NO NEED IN FIX FOR CURE53 ISSUE
                                        'stack_buffer_overflow__mark_alt_bonds_and_taut_groups'
                                        AS IT IS COVERED BY ALREADY DEFINED
                                        FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO
                                        */
#define FIX_GAF_2019_3 1
#define FIX_ONE_LINE_INCHI_INPUT_CONVERSION_ISSUE 1
#endif

#define ALLOW_EMPTY_INCHI_AS_INPUT    1 /* Allow "InChI=1//" and standard/beta analogues in input   */


#if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_USING_API) )
#define I2S_MODIFY_OUTPUT             1  /* 1=> Allow various InChI2InChI output types from cInChI */
#else
#define I2S_MODIFY_OUTPUT             0  /* 0=> Always */
#endif


#define FIX_NP_MINUS_BUG  1         /* 2010-03-11 DT Fix for bug reported by Timo Boehme						*/
                                    /* in normalization procedure for some structures containing N2(+) fragment	*/ 
                                    /* which may result in producing different InChI strings for the same		*/
                                    /* molecule, depending on original order of the atomic numbers				*/

/**************************/
/* additions to v1.00     */
/**************************/
#define FIX_ADJ_RAD                 0

#define SDF_OUTPUT_V2000            1  /* 1=>always output V2000 SDfile, 0=>only if needed */
#define SDF_OUTPUT_DT               1  /* 1=> all option -SdfAtomsDT to output D and T into SDfile */
#define CHECK_AROMBOND2ALT          1  /* 1=> check whether arom->alt bond conversion succeeded */

#define READ_INCHI_STRING           1  /* 1=> input InChI string and process it */
#if 0
#ifdef TARGET_LIB_FOR_WINCHI
#define READ_INCHI_STRING           0  /* 1=> input InChI string and process it */
#else
#define READ_INCHI_STRING           1  /* 1=> input InChI string and process it */
#endif
#endif

/****************************************************/
/* disabled extra external calls to InChI algorithm */
/****************************************************/
#define INCLUDE_NORMALIZATION_ENTRY_POINT  0

/**************************/
/* Normalization settings */
/**************************/

/* post version 1 features */
#define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
#define TAUT_15_NON_RING           1 /* 1,5 tautomerism with endpoints not in ring */

#define TAUT_PT_22_00              1 /* tautomerism rule PT_22_00 */
#define TAUT_PT_16_00              1 /* tautomerism rule PT_16_00 */
#define TAUT_PT_06_00              1 /* tautomerism rule PT_06_00 */
#define TAUT_PT_39_00              1 /* tautomerism rule PT_39_00 */
#define TAUT_PT_13_00              1 /* tautomerism rule PT_13_00 */
#define TAUT_PT_18_00              1 /* tautomerism rule PT_18_00 */  

#ifdef  BUILD_WITH_ENG_OPTIONS
#define UNDERIVATIZE               1 /* split to possible underivatized fragments */
#define RING2CHAIN                 1 /* open rings R-C(-OH)-O-R => R-C(=O) OH-R   */
#endif


#if( UNDERIVATIZE == 1 )
#define UNDERIVATIZE_REPORT        1 /* if SdfValue found, add to SdfValue a list of removed deriv. agents */

#define FIX_UNDERIV_TO_SDF       /* prevent bond normalization if underivatization result goes to SDF 2013-05-10 DCh */

#if 0
/*commented out 2019-08-20 as this switch caused several thousands changed InChI's for PubChem Substance (280M records) due to
[falsely] detected stereos in boranes, etc. e.g. 431803, 7890546,
*/
#define ALLOW_NO_CHARGE_ON_STEREO_CENTERS  /* do not require (+) on >N< for stereo, etc -- NIST output specific */
#endif

/*#define UNDERIV_SYLYL_ONLY */         /* Underiv: special case: recognize Sylyl derivatives only; typically commented out */

                                        /* derivatives selection begin */

                                        /* begin disabled derivatizations */
#ifdef NEVER
                                        /*#define UNDERIV_ACETATE_CnF2np1*/       /* 1r2c1-3 R-C(=O)-O---CnF2n+1 => R-C(=O)-OH, n=1..3: DERIV_BRIDGE_O - not a derivative */
                                                                                  /* Methyltion - 1 */
#define UNDERIV_ACETATE_Me            /* 1r1c3   R-C(=O)-O---Me (RCOO_Me) => R-C(=O)-OH: DERIV_BRIDGE_O */
                                                                                  /* Ethylation -1 */
#define UNDERIV_ACETATE_Et            /* 1r1c4   R-C(=O)-O---Et (RCOO_Et) => R-C(=O)-OH: DERIV_BRIDGE_O */
                                                                                  /* Propanoate - 3 */
#define UNDERIV_RN_AcEt               /* 2r1c4   R-N(-X)--C(=O)Et => R-N(-X)H: DERIV_BRIDGE_tN, X is not H */
#define UNDERIV_RNH_AcEt              /* 2r1c4   R-NH--C(=O)Et => R-NH2: DERIV_BRIDGE_NH */
#define UNDERIV_RO_COX_Et             /* 3r1c4   RO-C(=O)Et => ROH: create alcohols from acetates  DERIV_RO_COX */
#endif /* NEVER */
                                                                                  /* end disabled */

                                                                                  /* Acetate - 3 */
#define UNDERIV_RN_AcMe               /* 2r1c3   R-N(-X)--C(=O)Me => R-N(-X)H: DERIV_BRIDGE_tN, X is not H */
#define UNDERIV_RNH_AcMe              /* 2r1c3   R-NH--C(=O)Me => R-NH2: DERIV_BRIDGE_NH */
#define UNDERIV_RO_COX_Me             /* 3r1c1   RO-C(=O)Me => ROH: create alcohols from acetates  DERIV_RO_COX */
                                                                                  /* Benzoate - 1 */
#define UNDERIV_RO_COX_BENZOATES             /* 3r1c2 create alcohols from benzoates DERIV_RO_COX */

#define UNDERIV_RO_COX_PENTAFLOUROBENZOATES  /* 3r1c3 create alcohols from pentafluorobenzoates DERIV_RO_COX -C(=O)C6F5*/
#define UNDERIV_OOB_nButyl                   /* 4r2c1 DERIV_RING_O_OUTSIDE_PRECURSOR: 5 at, n-Butyl */
#define UNDERIV_X_OXIME_TBDMS                /* 5r2c3 DERIV_X_OXIME: >C=N--O-TBDMS */
#define UNDERIV_X_OXIME_TMS                  /* 5r2c2 DERIV_X_OXIME: >C=N--O-TMS */
#define UNDERIV_PYRROLIDIDES                 /* 7r1c1 DERIV_RING2_PRRLDD_OUTSIDE_PRECUR */

                                            /* derivatives selection end*/

#define UNDERIV_ADD_EXPLICIT_H               /* Underiv: add removed explict H after underivatization */
                                             /*#define UNDERIV_ADD_D_TO_PRECURSOR */      /* Uncomment to add Deuterium to the precursor struct -- for debugging only */
                                             /* search for other Underivatize settings in ichinorm.c, lines ~802-894 */
#endif /* UNDERIVATIZE == 1 */


/* post-2004-04-27 features */
#define HAL_ACID_H_XCHG            1 /* allow iso H exchange to HX (X=halogen) and H2Y (Y=halcogen) */
#define CANON_FIXH_TRANS           1 /* produce canonical fixed-H transposition */
#define STEREO_WEDGE_ONLY          1 /* 1=> only pointed ends stereo bonds define stereo; 0=> both ends */

/* current new (with respect to v1.12 Beta) preprocessing */
#define REMOVE_ION_PAIRS_EARLY     1 /* 1=> new preprocessing: step 1 before disconnecting metals in fix_odd_things() */
#define REMOVE_ION_PAIRS_DISC_STRU 1 /* 1=> new post-preprocessing: remove charhes after metal disconnection */
#define REMOVE_ION_PAIRS_FIX_BONDS 1 /* 1=> step2: set unchangeable bonds around removed ion pairs */
#define S_VI_O_PLUS_METAL_FIX_BOND 1 /* 1=> count double bond M-O(+)=S  as O=S in S(VI) ans S(VIII) fixing bonds */
#define N_V_STEREOBONDS            1 /* 1=> detect stereobonds incident to N(V); 0 => don't */
/* for testing */
#define REMOVE_ION_PAIRS_ORIG_STRU 0 /* 0=> normal mode (default)
                                      * 1=> testing mode only: remove ion pairs from the original structure
                                      *     to save the changes in the output Molfile (/OutputSDF) or AuxInfo
                                      *     NIP=No Ion Pairs
                                      */
/* salts treatment */
#define DISCONNECT_SALTS            1  /* 1=>disconnect metal atoms from salts, 0=>dont */ /* djb-rwth: default 1 */
#define TEST_REMOVE_S_ATOMS         1  /* 1=>default: after merging into one group test &
                                        *    remove unreachable,
                                        * 0=> old version: test only before merging into one t-group */
#define CHARGED_SALTS_ONLY          1  /* 1=>(default)do not test far salts tautomerism if
                                        *     no negative charge(s) present */
#define BNS_PROTECT_FROM_TAUT       1  /* 1=> do not allow testing of bonds to acetyl or nitro */
#define BNS_MARK_EDGE_2_DISCONNECT  1  /* 1=> mark edge as temp forbidden instead of disconnection */

#define REPLACE_ALT_WITH_TAUT       1  /* 1 => replace alt bonds with tautomeric bonds in case of standard t-groups */
#define MOVE_CHARGES                1  /* 1 => take moveable charges into account */
#define NEUTRALIZE_ENDPOINTS        1  /* 1 => before checking whether an H is moveable make 2 endpoints neutral */
                                       /*      implemented only if CHECK_TG_ALT_PATH = 0, defined in ichi_bns.c  */
#define FIX_H_CHECKING_TAUT         1  /* 1 => Fix moveable H or (-) before checking if taut. exchange is possible */
#define ALWAYS_ADD_TG_ON_THE_FLY    1  /* 1 => disables radical calcellation by taut-charge movement */
#define IGNORE_SINGLE_ENDPOINTS     1  /* 1 => see FindAccessibleEndPoints() in INChITaut.c */

/* recently added -- begin */
#define INCL_NON_SALT_CANDIDATATES   1  /* 1=> allow H and (-) migrate between "acidic" O and
                                         *     other possible endpoints */
#define SALT_WITH_PROTONS            1  /* 1=> (new new) include proton migrarion C-SH, =C-OH, NH+ */
#define OPPOSITE_CHARGE_IN_CGROUP    1  /* 1=> allow N(-) in (+) c-group, 0=> disallow */
#define MOVE_PPLUS_TO_REMOVE_PROTONS 0  /* 0=> default; 1=> (disabled) add P/P+ charge group during
                                         *     'hard' proton removal */
#define ADD_MOVEABLE_O_PLUS          1  /* 1=> allow charges on O(+) to move */
/* recently added -- end */

#define DISCONNECT_METALS           1  /* make main layer disconnected */ /* djb-rwth: default 1 */
#define RECONNECT_METALS            0  /* 1=> by default add reconnected layer in case of coord.
                                        *     compound disconnection */ /* djb-rwth: default 0 */
#define CHECK_METAL_VALENCE         0  /* 1=> disconnect only metals that have abnormal valence */
#define bREUSE_INCHI                1  /* 1=> do not recalulate INChI for components in reconnected
                                        *     structure that are same as in the connected one */
#define OUTPUT_CONNECTED_METAL_ONLY 0  /* 0=> default; 1 => (debug) create only reconnected or
                                        *     initial struct. output */
#define EMBED_REC_METALS_INCHI      1  /* 1=> (default) output Reconnected embedded in Disconnected INChI;
                                        * 0=> separate output */

#define bOUTPUT_ONE_STRUCT_TIME     1  /* 1 => output each structure time (non-release only) */



/* constants and array sizes */

#define INCHI_NUM            2    /* = array size; member indexes: */
#define INCHI_BAS            0    /* 0 => disconnected or normal */
#define INCHI_REC            1    /* 1 => reconnected */

#define TAUT_NUM            2    /* = array size; member indexes: */
#define TAUT_NON            0    /* 0 => normal structure */
#define TAUT_YES            1    /* 1 => tautomeric */
#define TAUT_INI            2    /* 2 => intermediate tautomeric structure */
#define ALT_TAUT(X)        ((X)>TAUT_YES? TAUT_YES : 1-(X))  /* was (1-(X)) */

/* INChI output modes */
#define OUT_N1              0    /* non-tautomeric only */
#define OUT_T1              1    /* tautomeric if present otherwise non-tautomeric */
#define OUT_NT              2    /* only non-taut representations of tautomeric */
#define OUT_TN              3    /* tautomeric if present otherwise non-tautomeric;
                                    separately output non-taut representations of tautomeric if present */
#define OUT_NN              4    /* only non-taut representations: non-taut else tautomeric */

/* OUT_TN = OUT_T1 + OUT_NT */

/* stereo */

#define NEW_STEREOCENTER_CHECK      1    /* 1 => add new stereocenter categories (see bCanInpAtomBeAStereoCenter(...)) */
#define MIN_SB_RING_SIZE            8    /* do not assume stereo bonds in rings containing 3..MIN_SB_RING_SIZE-1 atoms */

#define REMOVE_KNOWN_NONSTEREO      1 /* 1=> check in advance known stereo to remove parities from non-stereogenic elements */
#define REMOVE_CALC_NONSTEREO       1 /* 1=> check new stereo numberings to remove parities from non-stereogenic elements */
#define PROPAGATE_ILL_DEF_STEREO    1 /* 1=> if at least one of the pair of constitutionally identical (far) neighbors */
                                      /*     (of the tested atom) has ill-defined stereo parity and another has any */
                                      /*     stereo parity then set the parity of the tested atom to ill-defined value. */

#define ONLY_DOUBLE_BOND_STEREO     0  /* 1=> no alt bond stereo, no taut. bond attachment to stereo bond */
                                       /* 0=> allow other definitions (below) to be active */
#define ONE_BAD_SB_NEIGHBOR         1  /* 1 => allow 1 "bad" bond type neighbor to a stereobond atom. 2004-06-02 */

/* more stereo settings */
#define BREAK_ONE_MORE_SC_TIE       1   /* break one more tie when comparing possible stereocenter neighbors */
#define BREAK_ALSO_NEIGH_TIE        0   /* post 1.12Beta 2004-08-20: if fixed neighbor has equ neighbors, fix the one with smaller canon. rank */
#define BREAK_ALSO_NEIGH_TIE_ROTATE 1   /* post 1.12Beta 2004-09-02: break the second in 2nd psition; 1 works, 0 does not (example:MFCD01085607) */

#define STEREO_CENTER_BONDS_NORM   1   /* set length of the bonds around a stereocenter = 1 before getting the parity  */
#define STEREO_CENTER_BOND4_NORM   0   /* set length of the added bond around a stereocenter = 1 before getting the parity  */
#define NORMALIZE_INP_COORD        0   /* 0=>keep unchanged, 1 => make atom coordinates integer by normalizing to avg bond len 20 */

/* recent stereo */
#define STEREO_WEDGE_ONLY          1 /* 1=> only pointed ends stereo bonds define stereo; 0=> both ends 1.12Beta */
#define CHECK_C2v_S4_SYMM          0 /* post-1.12Beta 1=> check if a stereocenter has C2v or S4 symmetry; 0=>old mode */

#define EQL_H_NUM_TOGETHER          1 /* 1=> output 1-3,5H2 intead of 1-3H2,5H2 (CT_MODE_EQL_H_TOGETHER)  */
#define ABC_CT_NUM_CLOSURES         1 /* 1=> in coinnections compressed format output decimal number of closures instead of '-' */

/* temporary fix */
#define SINGLET_IS_TRIPLET          1 /* 'singlet' means two electrons make a lone pair instead of 2 bonds
                                         its effect on valence is same as the effect of a triplet */

/* defug: find structures where canonical partition is different from equitable */
#define FIND_CANON_NE_EQUITABLE     0  /* 0=>normal mode */
                                       /* 1=> extract (set EXTR_FLAGS = (EXTR_CANON_NE_EQUITABLE)*/
                                       /*     set cmd line options: /onlynonTAUT /: /UNCHARGEDACIDS:1 /DISCONSALT:0 /MOVEPOS:0 /DISCONMETAL:0 */

/* Debug: definitions for the extraction of the structures to the problem file */

/* definition of the flags for structure extraction to the
   problem file (for debugging and non-standard searching) */
#define EXTR_KNOWN_USED_TO_REMOVE_PARITY  0x000001
#define EXTR_CALC_USED_TO_REMOVE_PARITY   0x000002
#define EXTR_2EQL2CENTER_TO_REMOVE_PARITY 0x000004
#define EXTR_HAS_ATOM_WITH_DEFINED_PARITY 0x000008
#define EXTR_REMOVE_PARITY_WARNING        0x000010
#define EXTR_SALT_WAS_DISCONNECTED        0x000020
#define EXTR_SALT_PROTON_MOVED            0x000040
#define EXTR_SALT_PROTON_MOVE_ERR_WARN    0x000080
#define EXTR_METAL_WAS_DISCONNECTED       0x000100
#define EXTR_METAL_WAS_NOT_DISCONNECTED   0x000200
#define EXTR_NON_TRIVIAL_STEREO           0x000400 /* (Inv != Abs stereo) && (parities can't be obtained by inverting them) */
#define EXTR_UNUSUAL_VALENCES             0x000800
#define EXTR_HAS_METAL_ATOM               0x001000
#define EXTR_TEST_TAUT3_SALTS_DONE        0x002000 /* non-oxygen t-points used to discover tautomerism of merged t-groups */
#define EXTR_CANON_NE_EQUITABLE           0x004000 /* find structures where canonical partition is different from equitable */
#define EXTR_HAS_PROTON_PN                0x008000 /* has movable H+ attached to N or P */
#define EXTR_HAS_FEATURE                  0x010000 /* found a feature */
#define EXTR_TAUT_TREATMENT_CHARGES       0x020000 /* tautomeric treatment of charges */
#define EXTR_TRANSPOSITION_EXAMPLES       0x040000 /* extract structures that have different mobile-H and fixed-H orders */

/* define conditions of structure extraction to the problem file */
#define EXTR_MASK                        0 /*EXTR_TAUT_TREATMENT_CHARGES*/ /*(EXTR_HAS_FEATURE)*/ /*(EXTR_UNUSUAL_VALENCES | EXTR_HAS_METAL_ATOM)*/ /* 0 to disable */
#define EXTR_FLAGS                       0 /*EXTR_TAUT_TREATMENT_CHARGES*/ /*(EXTR_HAS_FEATURE)*/ /*(EXTR_HAS_PROTON_PN)*/ /*(EXTR_UNUSUAL_VALENCES)*/ /*(EXTR_CANON_NE_EQUITABLE)*/ /*(EXTR_TEST_TAUT3_SALTS_DONE)*/ /*(EXTR_HAS_METAL_ATOM)*/ /* (EXTR_NON_TRIVIAL_STEREO)*/ /*(EXTR_METAL_WAS_DISCONNECTED)*/ /* (EXTR_REMOVE_PARITY_WARNING)*/ /*(EXTR_HAS_ATOM_WITH_DEFINED_PARITY) */



/* added tautomeric structures */

#define TAUT_TROPOLONE_7            1  /* 1=> tautomeric 7-member rings ON */
#define TAUT_TROPOLONE_5            1  /* 1=> taut. similar to tropolone, 5-member ring */
#define TAUT_4PYRIDINOL_RINGS       1  /* 1=> OH-C5H4N rings tautomerism */
#define TAUT_PYRAZOLE_RINGS         1  /* 1=> tautomerizm in pyrazole rings */
/* limitation on tautomerism detection: */
#define TAUT_IGNORE_EQL_ENDPOINTS   0  /* 0=> even though 2 endpoints belong to same t-group check
                                              them to find more alt bonds (new)
                                          1=> ignore and do not check (old mode) */
#define TAUT_RINGS_ATTACH_CHAIN     1  /* 1=> allow only chain attachments to tautomeric endpoints */
                                       /*     (except pyrazole, where is no tautomeric attachment) */
                                       /* 0=> allow taut. attachments from same ring system. Default=1 */

#define FIND_RING_SYSTEMS           1  /* 1 => find and mark ring systems, blocks, cut-vertices */
                                       /* Needed for 5- and 6-member ring tautomers and in other places */

#define FIND_RINS_SYSTEMS_DISTANCES 0  /* 1 => find ring system and atom distance from terminal */
#define USE_DISTANCES_FOR_RANKING   0  /* 1 => rank ring systems according to distances from terminal */

#define DISPLAY_RING_SYSTEMS        0  /* 1 => for debug only; displays: */
                                       /* "block no"/"ring system no"/"cut-vertex (num. intersecting blocks-1)" */
                                       /* instead of ranks */
/* consistency */

#if ( bRELEASE_VERSION==1 && bOUTPUT_ONE_STRUCT_TIME==1)
#undef bOUTPUT_ONE_STRUCT_TIME
#define bOUTPUT_ONE_STRUCT_TIME 0
#endif

/* consistency: bRELEASE_VERSION==1 needs FIND_RING_SYSTEMS=1 */
#if ( bRELEASE_VERSION==1 && FIND_RING_SYSTEMS!=1 )
#ifdef FIND_RING_SYSTEMS
#undef FIND_RING_SYSTEMS
#endif
#define FIND_RING_SYSTEMS 1
#endif

/* consistency: FIND_RINS_SYSTEMS_DISTANCES needs FIND_RING_SYSTEMS  */
#if ( FIND_RING_SYSTEMS != 1 )

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
#undef  FIND_RINS_SYSTEMS_DISTANCES
#define FIND_RINS_SYSTEMS_DISTANCES 0
#endif

#endif

/* consistency: USE_DISTANCES_FOR_RANKING and DISPLAY_RING_SYSTEMS need FIND_RINS_SYSTEMS_DISTANCES */
#if ( FIND_RINS_SYSTEMS_DISTANCES != 1 )

#if ( USE_DISTANCES_FOR_RANKING == 1 )
#undef  USE_DISTANCES_FOR_RANKING
#define USE_DISTANCES_FOR_RANKING 0
#endif

#if ( DISPLAY_RING_SYSTEMS == 1 )
#undef  DISPLAY_RING_SYSTEMS
#define DISPLAY_RING_SYSTEMS 0
#endif

#endif


#if ( FIND_RING_SYSTEMS==1 && (TAUT_TROPOLONE_7==1 || TAUT_TROPOLONE_5==1 || TAUT_4PYRIDINOL_RINGS==1 || TAUT_PYRAZOLE_RINGS) )
#define TAUT_OTHER 1
#else
#define TAUT_OTHER 0
#endif

#define APPLY_IMPLICIT_H_DOWN_RULE 0   /* 1=> if 3 non-H atoms around stereocenter are in same plane */
                                       /*     then add "down" hydrogen to obtain sterecenter oparity */
                                       /* 0=> Implicit H stereo is unknown if all bonds to 3 non-H atoms */
                                       /*     are in XY plane */
#define ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS 1 /* 1=> consider bond in an alternating circuit stereogenic */
                                                 /*     even though it has adjacent tautomeric atom(s) */

#define IGNORE_TGROUP_WITHOUT_H   1    /* ignore tautomeric groups containing charges only */

#if ( DISCONNECT_SALTS == 1 )
#define REMOVE_TGROUP_CHARGE      0    /* 0: do not remove charge information from tautomeric groups */
#else
#define REMOVE_TGROUP_CHARGE      1    /* 1: remove charge information from tautomeric groups */
#endif

#if ( REMOVE_TGROUP_CHARGE == 1 )
#define INCHI_T_NUM_MOVABLE  1
#else
#define INCHI_T_NUM_MOVABLE  2
#endif

/******************************************/
/*   define canonicalization modes here   */
/******************************************/

#define USE_AUX_RANKING        1 /* 1=> get auxiliary ranking to accelerate canonicalization of H layers */
#define USE_AUX_RANKING_ALL    1 /* 1=> include all vertices in CellGetMinNode() selection 0=> only vertices with highest ranks */

#define USE_ISO_SORT_KEY_HFIXED  0  /* 0=> normal mode: merge isotopic taut H to isotopic atom sorting key in
                                           taut H-fixed canonicalization;
                                       1=> add one more "string" iso_sort_Hfixed to the canonicalization */

/************************
  questionable behavior
 ************************/
#define REL_RAC_STEREO_IGN_1_SC  0 /* 1=> drop from InChI sp3 stereo in components that have a single stereocenter */
                                   /* 0=> old-old mode (all such sp3 stereo is in the Identifier) */
/* internal definitions; see also REQ_MODE_BASIC etc in ichi.h */
#define CMODE_CT                 0x000001
#define CMODE_ISO                0x000002
#define CMODE_ISO_OUT            0x000004 /* obsolete ? */
#define CMODE_STEREO             0x000008
#define CMODE_ISO_STEREO         0x000010
#define CMODE_TAUT               0x000020
#define CMODE_NOEQ_STEREO        0x000040 /* 5-24-2002: do not use stereo equivalence to accelerate */
#define CMODE_REDNDNT_STEREO     0x000080 /* 6-11-2002: do not check for redundant stereo elements */
#define CMODE_NO_ALT_SBONDS      0x000100 /* 6-14-2002: do not assign stereo to alternating bonds */
/* new 10-10-2003 */
#define CMODE_RELATIVE_STEREO    0x000200    /* REL All Relative Stereo */
#define CMODE_RACEMIC_STEREO     0x000400    /* RAC All Racemic Stereo */
#define CMODE_SC_IGN_ALL_UU      0x000800    /* IAUSC Ignore stereocenters if All Undef/Unknown */
#define CMODE_SB_IGN_ALL_UU      0x001000    /* IAUSC Ignore stereobonds if All Undef/Unknown */
/* end of 10-10-2003 */

/* external definitions */
#define CANON_MODE_CT         (CMODE_CT)
#define CANON_MODE_TAUT       (CMODE_CT|CMODE_TAUT)
#define CANON_MODE_ISO        (CMODE_CT|CMODE_ISO|CMODE_ISO_OUT)
#define CANON_MODE_STEREO     (CMODE_CT|CMODE_STEREO)
#define CANON_MODE_ISO_STEREO (CMODE_CT|CMODE_ISO|CMODE_ISO_OUT|CMODE_ISO_STEREO)

#define CANON_MODE_MASK       0x00FF  /* used to determine canonicalization mode */

/*************************************************
 * from d_norm.c
 */

/* implemented definitions for CT_ATOMID */
#define CT_ATOMID_DONTINCLUDE   1
#define CT_ATOMID_IS_INITRANK   2
#define CT_ATOMID_IS_CURRANK    3

/***************************************
 * canonicalization settings  I
 ***************************************/

#define CANON_TAUTOMERS                1  /* 1=> process tautomers */
#define HYDROGENS_IN_INIT_RANKS        1  /* 1=> include num_H in initial ranking */

#define DOUBLE_BOND_NEIGH_LIST    0  /* 1 => include double bond neighbor in NeighList 2 times */
#define INCL_NON_6AROM            1  /* 1 => mark all arom. bonds; 0=>mark arom. bonds only in 6-member rings */

#define CT_SMALLEST             /* minimal CT */

#define CT_NEIGH_SMALLER        /* in CT, include neighbors with smaller ranks */

#define CT_ATOMID          CT_ATOMID_IS_CURRANK /*CT_ATOMID_DONTINCLUDE */

#define CT_NEIGH_INCREASE               /* in CT, neighbors ranks increase  */

#define USE_SYMMETRY_TO_ACCELERATE 1   /*1 => for fast CT canonicalization, to avoid full enumeration */

/* dependent definitions due to settings */

#ifdef CT_SMALLEST
#define CT_GREATER_THAN    >
#define CT_INITVALUE      ~0
#define BEST_PARITY        1  /* odd */
#define WORSE_PARITY       2
#else
#define CT_GREATER_THAN    <
#define CT_INITVALUE       0
#define BEST_PARITY        2  /* even */
#define WORSE_PARITY       1
#endif

#ifdef CT_NEIGH_SMALLER
#define CT_NEIGH_SMALLER_THAN <
#else
#define CT_NEIGH_SMALLER_THAN >
#endif

/* verify corectness of dependent settings */
#if !defined( CT_ATOMID )
#error  You have to #define CT_ATOMID
#else
#if ( defined( CT_ATOMID ) && CT_ATOMID==CT_ATOMID_DONTINCLUDE )
#error  CT_DELIMITER should be #defined if CT_ATOMID is not included
#endif
#endif

/***************************************
 * canonicalization settings  II
 ***************************************/
/* from extr_ct.h */
#define ALL_ALT_AS_AROMATIC         1  /* 1 => all altrnate bonds (even in cyclooctateraene) treat as aromatic */
                                       /*      and set DOUBLE_BOND_NEIGH_LIST = 0 */
#define ANY_ATOM_IN_ALT_CYCLE       1  /* 1=> accept any atom in alternating bond circuit, 0=>only some */

#define EXCL_ALL_AROM_BOND_PARITY   0  /* 1 => any arom atom cannot belong to stereo bond. */
                                       /*      This has presedence over ADD_6MEMB_AROM_BOND_PARITY=1 */
                                       /* 0 => include arom bonds parities according to */
                                       /*      ADD_6MEMB_AROM_BOND_PARITY definition */

#if ( EXCL_ALL_AROM_BOND_PARITY == 0 )
#define ADD_6MEMB_AROM_BOND_PARITY  1  /* 1 => all arom bonds are stereo bonds */
                                       /* 0 => only those arom bonds which do not belong to */
                                       /*      6-member arom rings are stereo bonds */
#else
#define ADD_6MEMB_AROM_BOND_PARITY  0  /* 0 => standard; 1 => meaningless: ignore parities of non-6-member ring alt. bonds */
#endif

#define MAX_NUM_STEREO_BONDS      3
#define MAX_NUM_STEREO_BOND_NEIGH 3
#define MIN_NUM_STEREO_BOND_NEIGH 2

#define MAX_NUM_STEREO_ATOM_NEIGH 4
#define STEREO_AT_MARK            8 /* > MAX_NUM_STEREO_BONDS */

#if ( ONLY_DOUBLE_BOND_STEREO == 1 )  /* { */

#ifdef ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS
#undef ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS
#define ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS 0
#endif

#ifdef EXCL_ALL_AROM_BOND_PARITY
#undef EXCL_ALL_AROM_BOND_PARITY
#define EXCL_ALL_AROM_BOND_PARITY 1
#endif

#ifdef ADD_6MEMB_AROM_BOND_PARITY
#undef ADD_6MEMB_AROM_BOND_PARITY
#define ADD_6MEMB_AROM_BOND_PARITY 0
#endif

#endif /* } ONLY_DOUBLE_BOND_STEREO */

/* dependent definitions due to settings */
#if ( ALL_ALT_AS_AROMATIC == 1 && DOUBLE_BOND_NEIGH_LIST != 0 )
#undef DOUBLE_BOND_NEIGH_LIST
#define DOUBLE_BOND_NEIGH_LIST 0
#endif


/*************************************
 * Drawing
 */

#define DRAW_AROM_TAUT  1              /* 1=> draw distinct aromatic & tautomer bonds, 0=> don't */

/******************************************************/
/*       C O M M O N     D E F I N I T I O N S        */
/******************************************************/


/* input bTautFlags flags */
#define TG_FLAG_TEST_TAUT__ATOMS         0x00000001   /* find regular tautomerism */
#define TG_FLAG_DISCONNECT_SALTS         0x00000002   /* DISCONNECT_SALTS disconnect */
#define TG_FLAG_TEST_TAUT__SALTS         0x00000004   /* DISCONNECT_SALTS if possible find long-range H/(-) taut. on =C-OH, >C=O    */
#define TG_FLAG_MOVE_POS_CHARGES         0x00000008   /* MOVE_CHARGES allow long-range movement of N(+), P(+) charges           */
#define TG_FLAG_TEST_TAUT2_SALTS         0x00000010   /* TEST_REMOVE_S_ATOMS multi-attachement long-range H/(-) taut. on =C-OH, >C=O   */
#define TG_FLAG_ALLOW_NO_NEGTV_O         0x00000020   /* CHARGED_SALTS_ONLY=0 (debug) find long-range H-only tautomerism on =C-OH, >C=O */
#define TG_FLAG_MERGE_TAUT_SALTS         0x00000040   /* DISCONNECT_SALTS merge all "salt"-t-groups and other =C-OH into one t-group */

#define TG_FLAG_ALL_TAUTOMERIC          (TG_FLAG_TEST_TAUT__ATOMS| \
                                         TG_FLAG_TEST_TAUT__SALTS| \
                                         TG_FLAG_TEST_TAUT2_SALTS| \
                                         TG_FLAG_MERGE_TAUT_SALTS)

#define TG_FLAG_DISCONNECT_COORD         0x00000080   /* find "coord. centers" and disconnect them */
#define TG_FLAG_RECONNECT_COORD          0x00000100   /* reconnect disconnected "coord. centers" */
#define TG_FLAG_CHECK_VALENCE_COORD      0x00000200   /* do not disconnect "coord. centers" with usual valence */
#define TG_FLAG_MOVE_HPLUS2NEUTR         0x00000400   /* move protons to neutralize */
#define TG_FLAG_VARIABLE_PROTONS         0x00000800   /* add/remove protons to neutralize */
#define TG_FLAG_HARD_ADD_REM_PROTONS     0x00001000   /* add/remove protons to neutralize in hard way */
#define TG_FLAG_POINTED_EDGE_STEREO      0x00002000   /* only pointed edge of stereo bond defines stereo */
#if ( FIX_ADJ_RAD == 1 )
#define TG_FLAG_FIX_ADJ_RADICALS         0x00004000   /* remove adjacent radical-doubletes, fix valence */
#endif
#define TG_FLAG_PHOSPHINE_STEREO         0x00008000   /* add phosphine sp3 stereo */
#define TG_FLAG_ARSINE_STEREO            0x00010000   /* add arsine sp3 stereo */
#define TG_FLAG_H_ALREADY_REMOVED        0x00020000   /* processing structure restored from InChI */
#define TG_FLAG_FIX_SP3_BUG              0x00040000   /* fix sp3 stereo bug: overlapping 2D stereo bond & coordinate scaling */

#define TG_FLAG_KETO_ENOL_TAUT           0x00080000   /* turn on keto-enol tautomerism detection */
#define TG_FLAG_1_5_TAUT                 0x00100000   /* turn on 1,5 tautomerism detection */

/* FB2 */
#define TG_FLAG_FIX_ISO_FIXEDH_BUG       0x00200000   /* fix bug found after v.102b (isotopic H representation)  */
#define TG_FLAG_FIX_TERM_H_CHRG_BUG      0x00400000   /* fix bug found after v.102b (moving H charge in 'remove_terminal_HDT') */

#define TG_FLAG_PT_22_00                 0x00800000
#define TG_FLAG_PT_16_00                 0x01000000
#define TG_FLAG_PT_06_00                 0x02000000
#define TG_FLAG_PT_39_00                 0x04000000
#define TG_FLAG_PT_13_00                 0x08000000
#define TG_FLAG_PT_18_00                 0x10000000

/* output bTautFlags flags */

#define TG_FLAG_MOVE_HPLUS2NEUTR_DONE    0x00000001   /* protons have been moved to neutralize */
#define TG_FLAG_TEST_TAUT__ATOMS_DONE    0x00000002
#define TG_FLAG_DISCONNECT_SALTS_DONE    0x00000004
#define TG_FLAG_TEST_TAUT__SALTS_DONE    0x00000008   /* multiple H tautomerism */
#define TG_FLAG_MOVE_POS_CHARGES_DONE    0x00000010
#define TG_FLAG_TEST_TAUT2_SALTS_DONE    0x00000020   /* merged t-groups */
#define TG_FLAG_ALLOW_NO_NEGTV_O_DONE    0x00000040
#define TG_FLAG_MERGE_TAUT_SALTS_DONE    0x00000080   /* added non-taut O to taut groups */

#define TG_FLAG_ALL_SALT_DONE          (TG_FLAG_TEST_TAUT__SALTS_DONE | \
                                        TG_FLAG_TEST_TAUT2_SALTS_DONE | \
                                        TG_FLAG_MERGE_TAUT_SALTS_DONE )

#define TG_FLAG_DISCONNECT_COORD_DONE    0x00000100   /* found and disconnected "coord. centers" */
#define TG_FLAG_CHECK_VALENCE_COORD_DONE 0x00000200   /* did not disconnect "coord. centers" with usual valence */
#define TG_FLAG_MOVE_CHARGE_COORD_DONE   0x00000400   /* changed charge of a disconnected ligand to fit its valence */
#define TG_FLAG_FIX_ODD_THINGS_DONE      0x00000800   /* fixed drawing ambiguities in fix_odd_things */
#define TG_FLAG_TEST_TAUT3_SALTS_DONE    0x00001000   /* merged t-groups + non-O taut atoms */
#define TG_FLAG_FOUND_SALT_CHARGES_DONE  0x00002000   /* not assigned: preprocessing detected possibility of salt-type tautomerism */
#define TG_FLAG_FOUND_ISOTOPIC_H_DONE    0x00004000   /* preprocessing detected isotopic H on "good" heteroatoms or isotopic H(+) */
#define TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE 0x00008000   /* preprocessing detected isotopic H on "good" heteroatoms or isotopic H(+) */
#if ( FIX_ADJ_RAD == 1 )
#define TG_FLAG_FIX_ADJ_RADICALS_DONE    0x00010000
#endif

#if ( READ_INCHI_STRING == 1 )
#define READ_INCHI_OUTPUT_INCHI          0x00000001
#define READ_INCHI_SPLIT_OUTPUT          0x00000002
#define READ_INCHI_KEEP_BALANCE_P        0x00000004
#define READ_INCHI_TO_STRUCTURE          0x00000008
#endif




/*********/
/*       */
/*  I/O  */
/*       */
/*********/

    typedef struct tagOutputString
    {
        char *pStr;
        int  nAllocatedLength;
        int  nUsedLength;
        int  nPtr;              /* expansion increment */
    } INCHI_IOS_STRING;


    /* INCHI_IOSTREAM.type values */
#define INCHI_IOS_TYPE_NONE 0
#define INCHI_IOS_TYPE_STRING 1
#define INCHI_IOS_TYPE_FILE 2

    typedef struct tagOutputStream
    {
        INCHI_IOS_STRING s;     /* output is directed either to resizable string buffer s   */
        FILE             *f;    /* or to the plain file:                                    */
        int              type;  /* dependent on type                                        */
    } INCHI_IOSTREAM;




    /***********/
    /*         */
    /*  DEBUG  */
    /*         */
    /***********/

#if ( defined(_WIN32) && defined(_DEBUG) && defined(_MSC_VER) /*&& !defined(COMPILE_ANSI_ONLY)*/ )
/* debug: memory leaks tracking */
#ifndef TARGET_LIB_FOR_WINCHI
#ifndef DO_NOT_TRACE_MEMORY_LEAKS
#define TRACE_MEMORY_LEAKS          1    /* 1=>trace, 0 => do not trace (Debug only) */
#else
#define TRACE_MEMORY_LEAKS          0
#endif
#else
#define TRACE_MEMORY_LEAKS          1    /* 1=>trace, **ALWAYS** =1 for TARGET_LIB_FOR_WINCHI */
#endif
#else /* not MSC and not Debug */
#define TRACE_MEMORY_LEAKS          0    /* 0: do not change */
#endif


/* memory leaks tracking */
#define INCHI_HEAPCHK          /* default: no explicit heap checking during the execution */

#if ( TRACE_MEMORY_LEAKS == 1 )
#ifdef _DEBUG

#define   inchi_malloc(s)         _malloc_dbg(s, _NORMAL_BLOCK, __FILE__, __LINE__)
#define   inchi_calloc(c, s)      _calloc_dbg(c, s, _NORMAL_BLOCK, __FILE__, __LINE__)
#define   inchi_realloc(d, s)     _realloc_dbg(d, s, _NORMAL_BLOCK, __FILE__, __LINE__)
#define   inchi_free(p)           _free_dbg(p, _NORMAL_BLOCK)

#ifdef TARGET_EXE_USING_API
/* INChI_MAIN specific */
#define e_inchi_malloc(a)   inchi_malloc(a)
#define e_inchi_calloc(a,b) inchi_calloc(a,b)
#define e_inchi_realloc(a,b) inchi_realloc(a,b)
#define e_inchi_free(a)     inchi_free(a)
#endif

/*#define  _CRTDBG_MAP_ALLOC*/  /* standard VC++ tool -- does not work with inchi_malloc(), etc */

#include <crtdbg.h>

/* to enable heap checking: #define CHECK_WIN32_VC_HEAP above #include "mode.h" in each source file or here */
#ifdef CHECK_WIN32_VC_HEAP
/* -- Confirms the integrity of the memory blocks allocated in the debug heap  -- */
#undef INCHI_HEAPCHK
#define INCHI_HEAPCHK \
do { \
    int tmp = _crtDbgFlag; \
    _crtDbgFlag |= _CRTDBG_ALLOC_MEM_DF; \
    _ASSERT( _CrtCheckMemory( ) ); \
    _crtDbgFlag = tmp; \
} while(0);

/*  -- less thorough than _CrtCheckMemory() check: check minimal consistency of the heap -- */
/*
#include <malloc.h>
#define INCHI_HEAPCHK \
do {\
   int heapstatus = _heapchk(); \
   _ASSERT( heapstatus != _HEAPBADBEGIN && heapstatus != _HEAPBADNODE && heapstatus != _HEAPBADPTR); \
} while(0);
*/
#endif

#else
#undef  TRACE_MEMORY_LEAKS
#define TRACE_MEMORY_LEAKS 0
#endif  /* _DEBUG */
#endif  /* TRACE_MEMORY_LEAKS */



/***********/
/*         */
/*  ALLOC  */
/*         */
/***********/

/* djb-rwth: fixing GH issue #89 */
#ifdef TARGET_EXE_USING_API
/* INChI_MAIN specific */
#ifndef inchi_malloc
#define inchi_malloc   e_inchi_malloc
#endif
#ifndef inchi_calloc
#define inchi_calloc   e_inchi_calloc
#endif
#ifndef inchi_realloc
#define inchi_realloc   e_inchi_realloc
#endif
#ifndef inchi_free
#define inchi_free     e_inchi_free
#endif

#ifndef e_inchi_malloc
#define e_inchi_malloc malloc
#endif
#ifndef e_inchi_calloc
#define e_inchi_calloc calloc
#endif
#ifndef e_inchi_realloc
#define e_inchi_realloc realloc
#endif
#ifndef e_inchi_free
#define e_inchi_free(X) do{ if(X) free(X); }while(0)
#endif

#else /* not TARGET_EXE_USING_API */

#ifndef inchi_malloc
#define inchi_malloc   malloc
#endif
#ifndef inchi_calloc
#define inchi_calloc   calloc
#endif
#ifndef inchi_realloc
#define inchi_realloc   realloc
#endif
#ifndef inchi_free
#define inchi_free(X)  do{ if(X) free(X); }while(0)
#endif

#endif /* TARGET_EXE_USING_API */

/* allocation/deallocation */
#define USE_ALLOCA 0

#if ( USE_ALLOCA == 1 )
#define qmalloc(X) _alloca(X)
#define qfree(X)   do{(X)=NULL;}while(0)
#else
#define qmalloc(X) inchi_malloc(X)
#define qfree(X)   do{if(X){inchi_free(X);(X)=NULL;}}while(0)
#endif

#if ( defined(_MSC_VER) && _MSC_VER >= 800 )
#define fast_alloc(X) _alloca(X)
#define fast_free(X)
#else
#define fast_alloc(X) inchi_malloc(X)
#define fast_free(X)  inchi_free(X)
#endif

#define qzfree(X)   do{if(X){inchi_free(X);(X)=NULL;}}while(0)

/* rellocation */
/* djb-rwth: avoiding memory leaks */
#define MYREALLOC2(PTRTYPE1, PTRTYPE2, PTR1, PTR2, LEN1, LEN2, ERR) \
    do { \
        if( (LEN1) <= (LEN2) ) {\
            PTRTYPE1 * newPTR1 = (PTRTYPE1 *)inchi_calloc( (LEN2)+1, sizeof(PTRTYPE1) );\
            PTRTYPE2 * newPTR2 = (PTRTYPE2 *)inchi_calloc( (LEN2)+1, sizeof(PTRTYPE2) );\
            if ( newPTR1 && newPTR2 ) { \
                if ( (PTR1) && (LEN1) > 0 ) \
                    (memcpy) ( newPTR1, (PTR1), (LEN1) * sizeof(PTRTYPE1) ); \
                if ( (PTR2) && (LEN1) > 0 ) \
                    (memcpy) ( newPTR2, (PTR2), (LEN1) * sizeof(PTRTYPE2) ); \
                if ( PTR1 ) \
                    inchi_free(PTR1);  \
                if ( PTR2 ) \
                    inchi_free(PTR2);  \
                (PTR1) = newPTR1; \
                (PTR2) = newPTR2; \
                (LEN1) = (LEN2);  \
                (ERR)  = 0; \
            } else {        \
                inchi_free(newPTR1); \
                inchi_free(newPTR2); \
                (ERR)  = 1; \
            }               \
        } else { (ERR) = 0; } \
    } while(0)


/* Solely v. 1.06 specific */

/* comment out the next line to enable polymeric debug */
#define DEBUG_POLYMERS 0
#ifndef DEBUG_POLYMERS
#define DEBUG_POLYMERS  2
#endif

#define POLYMERS_NO 0		/* ignore polymers										*/
#define POLYMERS_MODERN 1	/* v. 1.06+ way to treat polymers with Zz			*/
#define POLYMERS_LEGACY 2	/* v. 1.05 mode, no explicit Zz (internally they are here) */
#define POLYMERS_LEGACY_PLUS 3	/* v. 1.05 mode with an addition of that in all
                                   frame-shiftable-bistar-CRUs their backbone bonds 
                                   are reordered in descending seniority order. 
                                   Used as hidden 1st pass in 1.06 treatment		*/



#define STEREO_AT_ZZ 0
/* #define STEREO_AT_ZZ 1 */

/* Set to 1 and use in engineering mode, if necessary */
#define ALLOW_SUBSTRUCTURE_FILTERING 1


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /*  _MODE_H_ */
