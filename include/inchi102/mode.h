/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * January 10, 2009
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */




#ifndef __MODE_H__
#define __MODE_H__

#include <stdio.h>

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/*
    BUILD TARGET IS DETERMINED BY THE FOLLOWING DEFINES:

if (defined)        compile/build
-------------       -------------       

INCHI_LIBRARY       library for using INChI API described in inchi_api.h

INCHI_LINK_AS_DLL   INChI library as a Win32 DLL or to eliminate stricmp duplication

INCHI_MAIN          INCHI_MAIN.exe that calls INCHI_DLL.dll 

BUILD_CINCHI_WITH_INCHIKEY      
                    cInChI supporting InChIKey

INCHI_LIB           wInChI


*/


#define INCLUDE_V102_FEATURES 1
#ifdef INCLUDE_V102_FEATURES

#define BUILD_CINCHI_WITH_INCHIKEY 1
/* Comment the next line to use old (cInChI-1) generation defaults plus STDINCHI option */
#define STDINCHI_AS_DEFAULT 1
/*^^^ */
/* Comment the next line if not std-only package should be created */
#define ONLY_STDINCHI_STDKEY 1
#endif

/*^^^ */
#ifdef ONLY_STDINCHI_STDKEY
#define TREAT_UNKNOWN_STEREO_AS_UNDEFINED 1
#endif

/*^^^ */


/* uncomment to unconditionally force ANSI-89 C, no Win32 specific code */
#define INCHI_ANSI_ONLY

#ifdef INCHI_ANSI_ONLY
#define ADD_NON_ANSI_FUNCTIONS */ /* uncomment to add stricmp(), etc., see util.c */
#endif




/* uncomment to create a library for using INChI API described in inchi_api.h */
#define INCHI_LIBRARY

/* uncomment to use INChI library as a Win32 DLL or to eliminate stricmp duplication */
/* #define INCHI_LINK_AS_DLL */

/* Uncomment to compile INCHI_MAIN.exe that calls INCHI_DLL.dll */
/* #define INCHI_MAIN */

#if(!defined(_MSC_VER) || defined(INCHI_LIBRARY))  /* non-Microsoft GNU C, BCC, etc. compilers */
#ifndef INCHI_ANSI_ONLY
#define INCHI_ANSI_ONLY
#endif
#endif

/* #define INCHI_ALL_CPP */   /* uncomment to allow C++ compilation/linkage of functions prototyped in .h files */

#ifdef _MSC_VER
/*
========== disable MS VC++ 6.0 Level 4 compiler warnings: ==============
 C4706: assignment within conditional expression
 C4127: conditional expression is constant
 C4244: '=' : conversion from 'int ' to '???', possible loss of data
 C4701: local variable '???' may be used without having been initialized (removed)
 C4514: unreferenced inline/local function has been removed (C++)
 C4100: 'identifier' : unreferenced formal parameter
 C4786: 'identifier' : identifier was truncated to 'number' characters in the debug information
 C4996: 'identifier' was declared deprecated
========================================================================
*/
   #pragma warning( disable : 4706 4127 4514 4100 4786 4996 )
#endif

#ifndef SPECIAL_BUILD
#define SPECIAL_BUILD 0
#endif

#if( SPECIAL_BUILD == 1 )
#define SPECIAL_BUILD_STR  " (Modified Chemical Identifier - For Testing Only.)"
#endif

#ifdef  SPECIAL_BUILD_STR

#define SPECIAL_BUILD_STRING SPECIAL_BUILD_STR

#else

/*^^^ */
/* #define SPECIAL_BUILD_STRING ", Software version 1.01 release 07/21/2006" */
#define SPECIAL_BUILD_STRING ", Software version 1.02 release 01/10/2009"
/*^^^ */

#endif

/* Release */
#define bRELEASE_VERSION  1    /* 1=> release version; comment out to disable */
/* Debug */
#ifndef bRELEASE_VERSION 
#define bRELEASE_VERSION  0    /* 0=> debug version */
#endif

#define ACD_LABS_VERSION 0    /* always 0 */

#ifndef ADD_CMLPP
/* this allows ADD_CMLPP be #defined in a makefile */
#define ADD_CMLPP        0    /* 1 => add CMLPP input */
#endif

#if ( ADD_CMLPP == 1 )
#ifdef USE_CMLPPDLL
/* 1200 is VC++ 6.0 version, 1300 is VC++ .NET; USE_CMLPPDLL may be #defined in a makefile*/
#if ( defined(_WIN32) && defined(_MSC_VER) && _MSC_VER >= 1200 )
#define MSC_DELAY_LOAD_CMLPPDLL
#endif
#endif
#endif

/* display (non-canonical) c-groups, display orig at numbers */
#if( bRELEASE_VERSION == 1 )
#define DISPLAY_DEBUG_DATA_C_POINT 0  /* disabled release version for now */
#define DISPLAY_ORIG_AT_NUMBERS    1  /* 1 => in an uncanonicalized components display orig. atom numbers (default) */
#else
#define DISPLAY_DEBUG_DATA_C_POINT 1  /* debug: 1=>display (non-canonically numbered) c-groups, 0=>do not display */
#define DISPLAY_ORIG_AT_NUMBERS    1  /* 0 => in an uncanonicalized components display ordering atom numbers (debug) */
#endif

#if ( DISPLAY_DEBUG_DATA_C_POINT > 0 )
#define DISPLAY_DEBUG_DATA         DISPLAY_DEBUG_DATA_C_POINT
#endif

#define INCLUDE_V1_FIXES         1 /*SPECIAL_BUILD*/ /* prohibit all fixes to get v1 exact behavior */

#if( INCLUDE_V1_FIXES == 1 || SPECIAL_BUILD == 1 )

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

#else

#define FIX_ChCh_STEREO_CANON_BUG     0
#define ADD_ChCh_STEREO_CANON_CHK     0
#define FIX_ChCh_CONSTIT_CANON_BUG    0
#define FIX_EITHER_STEREO_IN_AUX_INFO 0
#define FIX_NORM_BUG_ADD_ION_PAIR     0
#define FIX_REM_PROTON_COUNT_BUG      0
#define FIX_READ_AUX_MEM_LEAK         0
#define FIX_READ_LONG_LINE_BUG        0
#define FIX_N_V_METAL_BONDS_GPF       0
#define BNS_RAD_SEARCH                0

#endif

/*******************************/
/* bug fixes in post-v1.00     */
/*******************************/
#if ( SPECIAL_BUILD == 1 )

#define FIX_ODD_THINGS_REM_Plus_BUG   1  /* 1 => fix bug: remove (1+) or (2+) charge from the central atom; */
                                         /* see Table 1, types 4 & 5 in InChI Tech Manual */
#define FIX_N_MINUS_NORN_BUG          1  /* to avoid numbering dependence of the normalization in case of */
                                         /* >N(+)=-=-NH(-) added neutralization before removing positive charge */
#define FIX_CANCEL_CHARGE_COUNT_BUG   1  /* 1 => fix bug in counting neutralization: only the last result is */
                                         /* remembered; 0 => keep as it was in version 1 */
#define FIX_2D_STEREO_BORDER_CASE     1  /* 1=> fix bug: 2D sp3 Stereo uncertainty in case of overlapping neighbors */
#define FIX_REM_ION_PAIRS_Si_BUG      1  /* 1 => detect Si same way as C in replacing ions with higher bond  */
                                         /* oreders and in fixing bonds. see Tables 2 and 5 in the InChI Tech */
                                         /* Manual 0=> keep as it was in version 1. */
#define FIX_STEREO_SCALING_BUG        1  /* 1=> fix bug that makes sp3 stereo undefined in case of bonds length > 20 */
#define FIX_EMPTY_LAYER_BUG           0  /* output "e" if the component has no stereo in the current layer */
                                         /* AND has stereo in the preceding layer */
#define FIX_EITHER_DB_AS_NONSTEREO    0  /* 1=> do not count ALT bond wuth Unknown stereo as stereogenic to close ALT-circuit */
/* new behavior as compared to v1 */
#define FIX_BOND23_IN_TAUT            0  /* 1=> when normalization results in bond order 2-3 try to interpret the bond as double */
/* not finished/tested */
#define FIX_TACN_POSSIBLE_BUG         0  /* fix (t-group)--(atom)--([-]c-group) behavior */
#define FIX_KEEP_H_ON_NH_ANION        0  /* in (de)protonation, disallow paths that convert -NH(-) into =N(-) */
#define FIX_AVOID_ADP                 0  /* avoid aggressive deprotonation; may make the results depend on numbering */
/* may change InChI */
#define FIX_NUM_TG                    1  /* in case of 7 or less atoms not enough BNS vertices allocated */
/* changes InChI for isothiocyanate */
#define FIX_CPOINT_BOND_CAP2          1  /* allow the leftmost bond in isothiocyanate N(-)=C=S become double */

#else

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



#define INCLUDE_V102_FEATURES 1

/* Ensure that all v. 1.02-final fixes and additions are in place */

#ifdef INCLUDE_V102_FEATURES

#ifndef ALLOW_STDINCHI_SHORTCUT 
#define ALLOW_STDINCHI_SHORTCUT       1  /*^^^ (Jan-Feb 2008) */
#endif
#ifndef FIX_OXOANION_DRAWING
#define FIX_OXOANION_DRAWING          1  /*^^^ post-1.02b */
#endif
#ifndef FIX_AMIDINIUM_DRAWING
#define FIX_AMIDINIUM_DRAWING         1  /*^^^ post-1.02b */
#endif
#ifndef FIX_TRANSPOSITION_CHARGE_BUG
#define FIX_TRANSPOSITION_CHARGE_BUG  1  /* (2008-01-02) fix bug that leads to missed charge in some cases when /o is present */
#endif

#ifndef FIX_TERM_H_CHRG_BUG
#define FIX_TERM_H_CHRG_BUG 1 /*^^^^^ (July 6, 2008 IPl) */ 
                              /* fix bug: in some cases (dependent on ordering numbers),
                                 moving a charge from terminal H to heavy atom resulted in
                                 neutralizing H but not adjusting charge of heavy atom */
#endif


#ifndef FIX_ISO_FIXEDH_BUG
#define FIX_ISO_FIXEDH_BUG            1 /* (2007-09-24) 1=> Fix bug: missing fixed-H iso segment in case of single removed D(+) */
#endif
#ifndef FIX_ISO_FIXEDH_BUG_READ
#define FIX_ISO_FIXEDH_BUG_READ       0 /* (2007-09-24) 1=> Accommodate this InChI bug in reading InChI */
#endif
#ifndef FIX_DALKE_BUGS
#define FIX_DALKE_BUGS                1  
#endif
#ifndef FIX_I2I_STEREOCONVERSION_BUG	
#define FIX_I2I_STEREOCONVERSION_BUG 1 /* (2008-03-06)   1=> Fix bug of i2i conversion SAbs-->(SRel||Srac) */
#endif
#ifndef FIX_I2I_STEREOCONVERSION_BUG2
#define FIX_I2I_STEREOCONVERSION_BUG2 1 /* (2008-04-02)   1=> Fix bug of i2i conversion (missed empty /t) */
#endif
#ifndef FIX_I2I_STEREOCONVERSION_BUG3
#define FIX_I2I_STEREOCONVERSION_BUG3 1 /* (2008-04-10)   1=> Fix bug of i2i conversion */
									    /* (missed repeating /s in FI after F for multi-component case) */
#endif


#endif


#if( !defined(INCHI_LIBRARY) && !defined(INCHI_MAIN) )
#define I2S_MODIFY_OUTPUT             1  /* 1=> Allow various InChI2InChI output types from cInChI */
#else
#define I2S_MODIFY_OUTPUT             0  /* 0=> Always */
#endif

#endif

/**************************/
/* additions to v1.00     */
/**************************/
#if ( SPECIAL_BUILD == 1 )
#define FIX_ADJ_RAD                 1  /* remove radicals from adjacent atoms, fix valences: -FIXRAD */
#define ADD_PHOSPHINE_STEREO        1  /* make >P- stereogenic atoms: -SPXYZ */
#define ADD_ARSINE_STEREO           1  /* make >As- stereogenic atoms: -SPXYZ */
#else
#define FIX_ADJ_RAD                 0
#define ADD_PHOSPHINE_STEREO        1
#define ADD_ARSINE_STEREO           1
#endif

#define SDF_OUTPUT_V2000            1  /* 1=>always output V2000 SDfile, 0=>only if needed */
#define SDF_OUTPUT_DT               1  /* 1=> all option -SdfAtomsDT to output D and T into SDfile */
#define CHECK_AROMBOND2ALT          1  /* 1=> check whether arom->alt bond conversion succeeded */

#ifdef INCHI_LIB
#define READ_INCHI_STRING           0  /* 1=> input InChI string and process it */
#else
#define READ_INCHI_STRING           1  /* 1=> input InChI string and process it */
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

/* for in-house use */
#define UNDERIVATIZE               0 /* split to possible underivatized fragments */
#define RING2CHAIN                 0 /* open rings R-C(-OH)-O-R => R-C(=O) OH-R   */

/* post-2004-04-27 features */
#define HAL_ACID_H_XCHG            1 /* allow iso H exchange to HX (X=halogen) and H2Y (Y=halcogen) */
#define CANON_FIXH_TRANS           1 /* produce canonical fixed-H transposition */
#define STEREO_WEDGE_ONLY          0 /* 1=> only pointed ends stereo bonds define stereo; 0=> both ends */

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
#define DISCONNECT_SALTS            1  /* 1=>disconnect metal atoms from salts, 0=>dont */
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

#define DISCONNECT_METALS           1  /* make main layer disconnected */
#define RECONNECT_METALS            0  /* 1=> by default add reconnected layer in case of coord.
                                        *     compound disconnection */  
#define CHECK_METAL_VALENCE         0  /* 1=> disconnect only metals that have abnormal valence */
#define bREUSE_INCHI                1  /* 1=> do not recalulate INChI for components in reconnected
                                        *     structure that are same as in the connected one */
#define OUTPUT_CONNECTED_METAL_ONLY 0  /* 0=> default; 1 => (debug) create only reconnected or
                                        *     initial struct. output */
#define EMBED_REC_METALS_INCHI      1  /* 1=> (default) output Reconnected embedded in Disconnected INChI;
                                        * 0=> separate output */

#define bOUTPUT_ONE_STRUCT_TIME     1  /* 1 => output each structure time (non-release only) */


/*#define INCHI_VERSION     "0.9Beta"  */
/*#define INCHI_VERSION     "0.91Beta" */  /* 10-10-2002: sent to Jonathan Goodman */
/*#define INCHI_VERSION     "0.92Beta" */  /* 11-15-2002: added Hill notation; sent to S.Heller & S.Stein */
/*#define INCHI_VERSION     "0.93Beta" */  /* 12-09-2002: Fixed isotopic canon. bug & chiralanes; sent to S.Heller & A. McNaught */
/*#define INCHI_VERSION     "0.931Beta" */  /* Non-BNS without salts released to PMR 01-2003 */
/*#define INCHI_VERSION     "0.932Beta" */      /* Released to CAS 04-01-2003:
                                            * - Improved taut. definitions as compared to 01-2003;
                                            * - fixed bug: non-isotopic components' stereo missing from isotopic stereo
                                            * - fixed bug: couldn't properly read Unix files (EOL = LF instead of CR/LF)
                                            *   (effective only for MS VC++, Borland and MinGW/GCC compiles that accept "rb" mode of fopen)
                                            *   DJGPP/GCC does not seem to need this fix.
                                            */
/*==== Release version ===*/
/*#define INCHI_VERSION     "0.94Beta" */       /* 02-27-2003: Balanced network search to find alt paths and non-stereo bonds;
                                                      Implemented salts disconnection; added (-) to taut groups */
/*#define INCHI_VERSION     "1.12Beta" */       /* 1.12: 07-06-2004: sort order: No H formula,..; Pointed end stereo ON, Aggressive (de)protonation OFF */
                                                /* 1.11: 05-19-2004: annotated plain text output, fixed bugs */
                                                /* 1.1:  04-08-2004: variable protonation version */
                                                /* 1.01: 12-23-2003 protected bonds, isotopic canonicalization in GetBaseCanonRanking() */
                                                /* 1.02: 01-26-2004 fixed new isotopic tgroup canon bug, molfile merge bug */

/*#define INCHI_VERSION       "1.0RC"*/             /* 02-07-2005 v1.0 Release Candidate */

#if ( SPECIAL_BUILD == 1 )

#define INCHI_VERSION       "1a"
#define INCHI_NAME          "MoChI"
#define INCHI_REC_NAME      "ReChI"

/* #define INCHI_VERSION       "1"     */
/* #define INCHI_NAME          "InChI" */

#else

#define INCHI_VERSION       "1"
#define INCHI_NAME          "InChI"
#define INCHI_REC_NAME      "ReChI"

#endif

#define INCHI_NAM_VER_DELIM "="

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

/* torture test */

#define TEST_RENUMB_ATOMS           0    /* 1 => heavy duty test by multiple renumbering of atoms */
#define TEST_RENUMB_NEIGH           1    /* 1 => randomly permutate neighbors */
#define TEST_RENUMB_SWITCH          0    /* 1 => display & output another (different) picture */
#define TEST_RENUMB_ATOMS_SAVE_LONGEST 0 /* 1 => save the component with largest processing time into the problem file */

#if ( defined(_WIN32) && defined(_DEBUG) && defined(_MSC_VER) /*&& !defined(INCHI_ANSI_ONLY)*/ )
/* debug: memory leaks tracking */
#ifndef INCHI_LIB

#ifndef DO_NOT_TRACE_MEMORY_LEAKS
#define TRACE_MEMORY_LEAKS          1    /* 1=>trace, 0 => do not trace (Debug only) */
#else
#define TRACE_MEMORY_LEAKS          0
#endif

#else

#define TRACE_MEMORY_LEAKS          1    /* 1=>trace, **ALWAYS** =1 for INCHI_LIB */

#endif

#else /* not MSC and not Debug */

#define TRACE_MEMORY_LEAKS          0    /* 0: do not change */

#endif

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
#define NORMALIZE_INP_COORD        0   /* 0=>keep unchanged, 1 => make atom coordinates integer values, avg bond len=20 */

/* recent stereo */
#define STEREO_WEDGE_ONLY          0 /* 1=> only pointed ends stereo bonds define stereo; 0=> both ends 1.12Beta */
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


#define ENTITY_REFS_IN_XML_MESSAGES 1 /* 1=> replace ' " < > & in error/warning messages with xml entity references */

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

#if( bRELEASE_VERSION==1 && bOUTPUT_ONE_STRUCT_TIME==1)
#undef bOUTPUT_ONE_STRUCT_TIME
#define bOUTPUT_ONE_STRUCT_TIME 0
#endif

/* consistency: bRELEASE_VERSION==1 needs FIND_RING_SYSTEMS=1 */
#if( bRELEASE_VERSION==1 && FIND_RING_SYSTEMS!=1 )
#ifdef FIND_RING_SYSTEMS
#undef FIND_RING_SYSTEMS
#endif
#define FIND_RING_SYSTEMS 1
#endif

/* consistency: FIND_RINS_SYSTEMS_DISTANCES needs FIND_RING_SYSTEMS  */
#if( FIND_RING_SYSTEMS != 1 )

#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
#undef  FIND_RINS_SYSTEMS_DISTANCES
#define FIND_RINS_SYSTEMS_DISTANCES 0
#endif

#endif

/* consistency: USE_DISTANCES_FOR_RANKING and DISPLAY_RING_SYSTEMS need FIND_RINS_SYSTEMS_DISTANCES */
#if( FIND_RINS_SYSTEMS_DISTANCES != 1 )

#if( USE_DISTANCES_FOR_RANKING == 1 )
#undef  USE_DISTANCES_FOR_RANKING
#define USE_DISTANCES_FOR_RANKING 0
#endif

#if( DISPLAY_RING_SYSTEMS == 1 )
#undef  DISPLAY_RING_SYSTEMS
#define DISPLAY_RING_SYSTEMS 0
#endif

#endif


#if( FIND_RING_SYSTEMS==1 && (TAUT_TROPOLONE_7==1 || TAUT_TROPOLONE_5==1 || TAUT_4PYRIDINOL_RINGS==1 || TAUT_PYRAZOLE_RINGS) )
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
#if( defined( CT_ATOMID ) && CT_ATOMID==CT_ATOMID_DONTINCLUDE )
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

#define CML_NUM_AT_IN_ATREF4      4
#define MAX_NUM_STEREO_BONDS      3
#define MAX_NUM_STEREO_BOND_NEIGH 3
#define MIN_NUM_STEREO_BOND_NEIGH 2

#define MAX_NUM_STEREO_ATOM_NEIGH 4
#define STEREO_AT_MARK            8 /* > MAX_NUM_STEREO_BONDS */

#if( ONLY_DOUBLE_BOND_STEREO == 1 )  /* { */

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
#if( ALL_ALT_AS_AROMATIC == 1 && DOUBLE_BOND_NEIGH_LIST != 0 )
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
#if( FIX_ADJ_RAD == 1 )
#define TG_FLAG_FIX_ADJ_RADICALS         0x00004000   /* remove adjacent radical-doubletes, fix valence */
#endif
#define TG_FLAG_PHOSPHINE_STEREO         0x00008000   /* add phosphine sp3 stereo */
#define TG_FLAG_ARSINE_STEREO            0x00010000   /* add arsine sp3 stereo */
#define TG_FLAG_H_ALREADY_REMOVED        0x00020000   /* processing structure restored from InChI */
#define TG_FLAG_FIX_SP3_BUG              0x00040000   /* fix sp3 stereo bug: overlapping 2D stereo bond & coordinate scaling */

#define TG_FLAG_KETO_ENOL_TAUT           0x00080000   /* turn on keto-enol tautomerism detection */
#define TG_FLAG_1_5_TAUT                 0x00100000   /* turn on 1,5 tautomerism detection */
               
/*^^^ FB2 */
#define TG_FLAG_FIX_ISO_FIXEDH_BUG       0x00200000   /* fix bug found after v.102b (isotopic H representation)  */
#define TG_FLAG_FIX_TERM_H_CHRG_BUG      0x00400000   /* fix bug found after v.102b (moving H charge in 'remove_terminal_HDT') */
                        
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
#if( FIX_ADJ_RAD == 1 )
#define TG_FLAG_FIX_ADJ_RADICALS_DONE    0x00010000
#endif

#if( READ_INCHI_STRING == 1 )
#define READ_INCHI_OUTPUT_INCHI          0x00000001
#define READ_INCHI_SPLIT_OUTPUT          0x00000002
#define READ_INCHI_KEEP_BALANCE_P        0x00000004
#define READ_INCHI_TO_STRUCTURE          0x00000008
#endif


#ifdef _WIN32

#define   INCHI_OPTION_PREFX  '/'
#define   INCHI_PATH_DELIM    '\\'

#else

#define   INCHI_OPTION_PREFX  '-'
#define   INCHI_PATH_DELIM    '/'

#endif

#define   INCHI_ALT_OPT_PREFIX  '-'
#define   INCHI_ACD_LABS_PREFIX '-'

typedef struct tagOutputString {
    char *pStr;
    int  nAllocatedLength;
    int  nUsedLength;
    int  nPtr;
} INCHI_OUTPUT;


/*^^^ Post-1.02b */
typedef struct tagOutputStream 
{
    /* output is directed either to resizable string buffer: */
    INCHI_OUTPUT s;
    /* or to the plain file: */
    FILE* f;
    int type;
} INCHI_IOSTREAM;
/* INCHI_IOSTREAM.type values */
#define INCHI_IOSTREAM_NONE 0
#define INCHI_IOSTREAM_STRING 1
#define INCHI_IOSTREAM_FILE 2


/* memory leaks tracking */
#define INCHI_HEAPCHK          /* default: no explicit heap checking during the execution */

#if( TRACE_MEMORY_LEAKS == 1 )
#ifdef _DEBUG

#define   inchi_malloc(s)         _malloc_dbg(s, _NORMAL_BLOCK, __FILE__, __LINE__)
#define   inchi_calloc(c, s)      _calloc_dbg(c, s, _NORMAL_BLOCK, __FILE__, __LINE__)
#define   inchi_free(p)           _free_dbg(p, _NORMAL_BLOCK)

#ifdef INCHI_MAIN
/* INChI_MAIN specific */
#define e_inchi_malloc(a)   inchi_malloc(a)
#define e_inchi_calloc(a,b) inchi_calloc(a,b)
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

#ifdef INCHI_MAIN
/* INChI_MAIN specific */
#ifndef inchi_malloc
#define inchi_malloc   e_inchi_malloc
#endif
#ifndef inchi_calloc
#define inchi_calloc   e_inchi_calloc
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
#ifndef e_inchi_free
#define e_inchi_free(X) do{ if(X) free(X); }while(0)
#endif

#else /* not INCHI_MAIN */

#ifndef inchi_malloc
#define inchi_malloc   malloc
#endif
#ifndef inchi_calloc
#define inchi_calloc   calloc
#endif
#ifndef inchi_free
#define inchi_free(X)  do{ if(X) free(X); }while(0)
#endif

#endif /* INCHI_MAIN */

/* allocation/deallocation */
#define USE_ALLOCA 0

#if( USE_ALLOCA == 1 )
#define qmalloc(X) _alloca(X)
#define qfree(X)   do{(X)=NULL;}while(0)
#else
#define qmalloc(X) inchi_malloc(X)
#define qfree(X)   do{if(X){inchi_free(X);(X)=NULL;}}while(0)
#endif

#if( defined(_MSC_VER) && _MSC_VER >= 800 )
#define fast_alloc(X) _alloca(X)
#define fast_free(X)
#else
#define fast_alloc(X) inchi_malloc(X)
#define fast_free(X)  inchi_free(X)
#endif

#define qzfree(X)   do{if(X){inchi_free(X);(X)=NULL;}}while(0)

/* rellocation */

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
                (ERR)  = 1; \
            }               \
        } else { (ERR) = 0; } \
    } while(0)

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __MODE_H__ */
