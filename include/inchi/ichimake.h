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


#ifndef _ICHIMAKE_H_
#define _ICHIMAKE_H_

#include "ichisize.h"
#include "ichi.h"
#include "extr_ct.h"
#include "ichicant.h"

/***********************************************************************/
/* replace all ' ' delimiters with ',' */
#define ITEM_DELIMETER    ","
#define EXTRA_SPACE       ""
#define COMMA_EXTRA_SPACE ","
#define LEN_EXTRA_SPACE   0

/**********************************************************************************************/
/*  nCtMode for output INChI */

#define CT_MODE_NO_ORPHANS        1 /* no orphans, that CT should have only atoms with neighbors */
#define CT_MODE_ABC_NUMBERS       2
#define CT_MODE_ATOM_COUNTS       4
#define CT_MODE_PREDECESSORS      8
#define CT_MODE_EQL_H_TOGETHER   16
#define CT_MODE_ABC_NUM_CLOSURES 32 /* in CT_MODE_ABC_NUMBERS output AB1AC2AB instead of AB-AC-A-B */


/*************** Macros for retrieving requested INChI and INChI_Aux *****************************/
/* S->pINChI[TAUT_YES] has info: */
#define HAS_T(S)  (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms)
/* S->pINChI[TAUT_NON] has info: */
#define HAS_N(S)  (S->pINChI[TAUT_NON] && S->pINChI[TAUT_NON]->nNumberOfAtoms)

/* S->pINChI[TAUT_YES] has tautomeric info: */
#define HAS_TT(S) (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms && S->pINChI[TAUT_YES]->lenTautomer>0)
/* S->pINChI[TAUT_YES] has non-taitomeric info: */
#define HAS_TN(S) (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms && !S->pINChI[TAUT_YES]->lenTautomer)
/* S->pINChI[TAUT_NON] has non-tautomeric info: */
#define HAS_NN(S) (S->pINChI[TAUT_NON] && S->pINChI[TAUT_NON]->nNumberOfAtoms && !S->pINChI[TAUT_NON]->lenTautomer)
#define GET_II(M,S) ((M==OUT_N1)?              (HAS_TN(S)? TAUT_YES : HAS_NN(S)? TAUT_NON : -1): \
                     (M==OUT_T1 || M==OUT_TN)? (HAS_T(S) ? TAUT_YES : HAS_N(S) ? TAUT_NON : -1): \
                     (M==OUT_NN)?              (HAS_NN(S)? TAUT_NON : HAS_TN(S)? TAUT_YES : -1): \
                     (M==OUT_NT)?              ((HAS_TT(S) && HAS_NN(S))       ? TAUT_NON : -1) : -1)

/*********************************/
/* Equivalence flags definitions */
/*********************************/

/* What is equal (ii = INChI_ITEM) */
#define iiSTEREO      0x0001  /* stereo (sp2 or sp3)   */
#define iiSTEREO_INV  0x0002  /* inverted stereo (sp3) */
#define iiNUMB        0x0004  /* numbering or inverted stereo numbering */
#define iiEQU         0x0008  /* equivalence info */
/* derived:
  (iiSTEREO_INV | iiNUMB) = numbering of inverted stereo
*/

/* Additional info to what is equal (INCHI_ITEM_TYPE = iit) */
#define iitISO       0x0010  /* Isotopic */
#define iitNONTAUT   0x0020  /* Non-tautomeric */
/* derived:
  (iitISO | iitNONTAUT) = isotopic non-tautomeric
*/

/* Where is the equivalent item located (INChI_ITEM_EQUAL_TO = iiEq2) */
#define iiEq2NONTAUT 0x0040  /* non-tautomeric */
#define iiEq2ISO     0x0080  /* isotopic */
#define iiEq2INV     0x0100  /* equal to inverted (stereo sp3) or to numbering of inverted stereo */

#define iiEmpty      0x0200  /* item is empty while in the preceding layer the item is not empty */

/*********************** Printing strings external declarations *******************************/

extern const char sCompDelim[];

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/**********************************************************************************************/
int CompareTautNonIsoPartOfINChI( const INChI *i1,
                                      const INChI *i2 );
const char *EquString( int EquVal );
    int FillOutINChI( INChI *pINChI,
                      INChI_Aux *pINChI_Aux,
                      int num_atoms,
                      int num_at_tg,
                      int num_removed_H,
                      sp_ATOM *at,
                      inp_ATOM *norm_at,
                      CANON_STAT *pCS,
                      CANON_GLOBALS *pCG,
                      int bTautomeric,
                      INCHI_MODE nUserMode,
                      char *pStrErrStruct,
                      int bNoWarnings );
    int MakeHillFormulaString( char *szHillFormula,
                               INCHI_IOS_STRING *strbuf,
                               int *bOverflow );
    int bHasOrigInfo( ORIG_INFO *OrigInfo,
                      int num_atoms );
    int EqlOrigInfo( INChI_Aux *a1,
                     INChI_Aux *a2 );
    int MakeAbcNumber( char *szString,
                       int nStringLen,
                       const char *szLeadingDelim,
                       int nValue );
    int MakeDecNumber( char *szString,
                       int nStringLen,
                       const char *szLeadingDelim,
                       int nValue );
    int MakeCtStringNew( CANON_GLOBALS *pCG,
                         AT_NUMB *LinearCT,
                         int nLenCT,
                         int bAddDelim,
                         S_CHAR *nNum_H,
                         int num_atoms,
                         INCHI_IOS_STRING *strbuf,
                         int nCtMode,
                         int *bOverflow );
    int MakeCtStringOld( AT_NUMB *LinearCT,
                         int nLenCT,
                         int bAddDelim,
                         INCHI_IOS_STRING *strbuf,
                         int nCtMode,
                         int *bOverflow );
    int MakeCtString( CANON_GLOBALS *pCG,
                      AT_NUMB *LinearCT,
                      int nLenCT,
                      int bAddDelim,
                      S_CHAR *nNum_H,    /* not used here */
                      int num_atoms,    /* not used here */
                      INCHI_IOS_STRING *strbuf,
                      int nCtMode,
                      int *bOverflow );
    int MakeTautString( AT_NUMB *LinearCT,
                        int nLenCT,
                        int bAddDelim,
                        INCHI_IOS_STRING *strbuf,
                        int nCtMode,
                        int *bOverflow );
    int MakeEquString( AT_NUMB *LinearCT,
                       int nLenCT,
                       int bAddDelim,
                       INCHI_IOS_STRING *strbuf,
                       int nCtMode,
                       int *bOverflow );
    int MakeIsoAtomString( INChI_IsotopicAtom   *IsotopicAtom,
                           int nNumberOfIsotopicAtoms,
                           INCHI_IOS_STRING *strbuf,
                           int nCtMode,
                           int *bOverflow );
    int MakeIsoTautString( INChI_IsotopicTGroup   *IsotopicTGroup,
                           int nNumberOfIsotopicTGroups,
                           INCHI_IOS_STRING *strbuf,
                           int nCtMode,
                           int *bOverflow );
    int MakeIsoHString( int num_iso_H[],
                        INCHI_IOS_STRING *strbuf,
                        int nCtMode, int *bOverflow );
    int MakeStereoString( AT_NUMB *at1,
                          AT_NUMB *at2,
                          S_CHAR *parity,
                          int bAddDelim,
                          int nLenCT,
                          INCHI_IOS_STRING *buf,
                          int nCtMode,
                          int *bOverflow );
    int MakeCRVString( ORIG_INFO *OrigInfo,
                       int nLenCT,
                       int bAddDelim,
                       INCHI_IOS_STRING *strbuf,
                       int nCtMode,
                       int *bOverflow );
    int MakeMult( int mult,
                  const char *szTailingDelim,
                  INCHI_IOS_STRING *buf,
                  int nCtMode,
                  int *bOverflow );
    int MakeDelim( const char *szTailingDelim,
                   INCHI_IOS_STRING *buf,
                   int *bOverflow );
    int MakeEqStr( const char *szTailingDelim,
                   int mult,
                   INCHI_IOS_STRING *buf,
                   int *bOverflow );
    int MakeHString( int bAddDelim,
                     S_CHAR *LinearCT,
                     int nLenCT,
                     INCHI_IOS_STRING *buf,
                     int nCtMode,
                     int *bOverflow );
    AT_NUMB  *GetDfsOrder4CT( CANON_GLOBALS *pCG,
                              AT_NUMB *LinearCT,
                              int nLenCT,
                              S_CHAR *nNum_H,
                              int num_atoms,
                              int nCtMode );
    int str_HillFormula( INCHI_SORT *pINChISort,
                         INCHI_IOS_STRING *strbuf,
                         int *bOverflow,
                         int bOutType,
                         int num_components,
                         int bUseMulipliers );
    int str_HillFormula2( INCHI_SORT *pINChISort    /* non-taut */,
                          INCHI_SORT *pINChISort2    /* taut */,
                          INCHI_IOS_STRING *strbuf,
                          int *bOverflow,
                          int bOutType,
                          int num_components,
                          int bUseMulipliers );
    int str_Connections( CANON_GLOBALS *pCG,
                         INCHI_SORT *pINChISort,
                         INCHI_IOS_STRING *strbuf,
                         int *bOverflow,
                         int bOutType,
                         int ATOM_MODE,
                         int num_components,
                         int bUseMulipliers );
    int str_H_atoms( INCHI_SORT *pINChISort,
                     INCHI_IOS_STRING *strbuf,
                     int *bOverflow,
                     int bOutType,
                     int ATOM_MODE,
                     int TAUT_MODE,
                     int num_components,
                     int bUseMulipliers );
    int str_Charge2( INCHI_SORT *pINChISort,
                     INCHI_SORT *pINChISort2,
                     INCHI_IOS_STRING *strbuf,
                     int *bOverflow,
                     int bOutType,
                     int num_components,
                     int bSecondNonTautPass,
                     int bOmitRepetitions,
                     int bUseMulipliers );
    int str_Sp2( INCHI_SORT *pINChISort,
                 INCHI_SORT *pINChISort2,
                 INCHI_IOS_STRING *strbuf,
                 int *bOverflow,
                 int bOutType,
                 int TAUT_MODE,
                 int num_components,
                 int bSecondNonTautPass,
                 int bOmitRepetitions,
                 int bUseMulipliers );
    int str_IsoSp2( INCHI_SORT *pINChISort,
                    INCHI_SORT *pINChISort2,
                    INCHI_IOS_STRING *strbuf,
                    int *bOverflow,
                    int bOutType, int TAUT_MODE,
                    int num_components,
                    int bSecondNonTautPass,
                    int bOmitRepetitions,
                    int bUseMulipliers );
    int str_Sp3( INCHI_SORT *pINChISort,
                 INCHI_SORT *pINChISort2,
                 INCHI_IOS_STRING *strbuf,
                 int *bOverflow,
                 int bOutType,
                 int TAUT_MODE,
                 int num_components,
                 int bRelRac,
                 int bSecondNonTautPass,
                 int bOmitRepetitions,
                 int bUseMulipliers );
    int str_IsoSp3( INCHI_SORT *pINChISort,
                    INCHI_SORT *pINChISort2,
                    INCHI_IOS_STRING *strbuf,
                    int *bOverflow,
                    int bOutType,
                    int TAUT_MODE,
                    int num_components,
                    int bRelRac,
                    int bSecondNonTautPass,
                    int bOmitRepetitions,
                    int bUseMulipliers );
    int str_StereoAbsInv( INCHI_SORT *pINChISort,
                          INCHI_IOS_STRING *strbuf,
                          int *bOverflow,
                          int bOutType,
                          int num_components );
    int str_IsoStereoAbsInv( INCHI_SORT *pINChISort,
                             INCHI_IOS_STRING *strbuf,
                             int *bOverflow,
                             int bOutType,
                             int num_components );
    int str_IsoAtoms( INCHI_SORT *pINChISort,
                      INCHI_SORT *pINChISort2,
                      INCHI_IOS_STRING *strbuf,
                      int *bOverflow,
                      int bOutType,
                      int TAUT_MODE,
                      int num_components,
                      int bAbcNumbers,
                      int bSecondNonTautPass,
                      int bOmitRepetitions,
                      int bUseMulipliers );
    int str_FixedH_atoms( INCHI_SORT *pINChISort,
                          INCHI_IOS_STRING *strbuf,
                          int *bOverflow,
                          int bOutType,
                          int ATOM_MODE,
                          int num_components,
                          int bUseMulipliers );
    int str_AuxNumb( CANON_GLOBALS *pCG,
                     INCHI_SORT *pINChISort,
                     INCHI_SORT *pINChISort2,
                     INCHI_IOS_STRING *strbuf,
                     int *bOverflow,
                     int bOutType,
                     int TAUT_MODE,
                     int num_components,
                     int bSecondNonTautPass,
                     int bOmitRepetitions );
    int str_AuxEqu( INCHI_SORT *pINChISort,
                    INCHI_SORT *pINChISort2,
                    INCHI_IOS_STRING *strbuf,
                    int *bOverflow,
                    int bOutType,
                    int TAUT_MODE,
                    int num_components,
                    int bSecondNonTautPass,
                    int bOmitRepetitions,
                    int bUseMulipliers );
    int str_AuxTgroupEqu( INCHI_SORT *pINChISort,
                          INCHI_IOS_STRING *strbuf,
                          int *bOverflow,
                          int bOutType,
                          int TAUT_MODE,
                          int num_components,
                          int bUseMulipliers );
    int str_AuxIsoTgroupEqu( INCHI_SORT *pINChISort,
                             INCHI_IOS_STRING *strbuf,
                             int *bOverflow,
                             int bOutType,
                             int TAUT_MODE,
                             int num_components,
                             int bOmitRepetitions,
                             int bUseMulipliers );
    int str_AuxInvSp3( INCHI_SORT *pINChISort,
                       INCHI_SORT *pINChISort2,
                       INCHI_IOS_STRING *strbuf,
                       int *bOverflow,
                       int bOutType,
                       int TAUT_MODE,
                       int num_components,
                       int bSecondNonTautPass,
                       int bOmitRepetitions,
                       int bUseMulipliers );
    int str_AuxInvSp3Numb( CANON_GLOBALS *pCG,
                           INCHI_SORT *pINChISort,
                           INCHI_SORT *pINChISort2,
                           INCHI_IOS_STRING *strbuf,
                           int *bOverflow,
                           int bOutType,
                           int TAUT_MODE,
                           int num_components,
                           int bSecondNonTautPass,
                           int bOmitRepetitions );
    int str_AuxIsoNumb( CANON_GLOBALS *pCG,
                        INCHI_SORT *pINChISort,
                        INCHI_SORT *pINChISort2,
                        INCHI_IOS_STRING *strbuf,
                        int *bOverflow,
                        int bOutType,
                        int TAUT_MODE,
                        int num_components,
                        int bSecondNonTautPass,
                        int bOmitRepetitions );
    int str_AuxIsoEqu( INCHI_SORT *pINChISort,
                       INCHI_SORT *pINChISort2,
                       INCHI_IOS_STRING *strbuf,
                       int *bOverflow,
                       int bOutType,
                       int TAUT_MODE,
                       int num_components,
                       int bSecondNonTautPass,
                       int bOmitRepetitions,
                       int bUseMulipliers );
    int str_AuxInvIsoSp3( INCHI_SORT *pINChISort,
                          INCHI_SORT *pINChISort2,
                          INCHI_IOS_STRING *strbuf,
                          int *bOverflow, int bOutType,
                          int TAUT_MODE,
                          int num_components,
                          int bSecondNonTautPass,
                          int bOmitRepetitions,
                          int bUseMulipliers );
    int str_AuxInvIsoSp3Numb( CANON_GLOBALS *pCG,
                              INCHI_SORT *pINChISort,
                              INCHI_SORT *pINChISort2,
                              INCHI_IOS_STRING *strbuf,
                              int *bOverflow,
                              int bOutType,
                              int TAUT_MODE,
                              int num_components,
                              int bSecondNonTautPass,
                              int bOmitRepetitions );
    int str_AuxChargeRadVal( INCHI_SORT *pINChISort,
                             INCHI_IOS_STRING *strbuf,
                             int *bOverflow,
                             int bOutType,
                             int TAUT_MODE,
                             int num_components,
                             int bUseMulipliers );
    int bin_AuxTautTrans( INCHI_SORT *pINChISort,
                          INCHI_SORT *pINChISort2,
                          AT_NUMB **pTrans_n,
                          AT_NUMB **pTrans_s,
                          int bOutType,
                          int num_components );
int str_AuxTautTrans( CANON_GLOBALS *pCG,
                      AT_NUMB *nTrans_n,
                      AT_NUMB *nTrans_s,
                      INCHI_IOS_STRING *strbuf,
                      int *bOverflow,
                      int TAUT_MODE,
                      int num_components );

int MergeZzInHillFormula(INCHI_IOS_STRING *strbuf);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



#endif    /* _ICHIMAKE_H_ */
