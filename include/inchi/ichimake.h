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


#ifndef __INCHIMAKE_H__
#define __INCHIMAKE_H__

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

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/**********************************************************************************************/
int FillOutINChI( INChI *pINChI, INChI_Aux *pINChI_Aux,
                 int num_atoms, int num_at_tg, int num_removed_H,
                 sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, int bTautomeric,
                 INCHI_MODE nUserMode, char *pStrErrStruct );



int MakeHillFormulaString( char *szHillFormula, char *szLinearCT, int nLen_szLinearCT, int *bOverflow);

int bHasOrigInfo( ORIG_INFO *OrigInfo, int num_atoms );
int EqlOrigInfo( INChI_Aux *a1, INChI_Aux *a2 );

int MakeAbcNumber( char *szString, int nStringLen, const char *szLeadingDelim, int nValue );
int MakeDecNumber( char *szString, int nStringLen, const char *szLeadingDelim, int nValue );
int MakeCtStringNew( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  S_CHAR *nNum_H, int num_atoms,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeCtStringOld( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeCtString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  S_CHAR *nNum_H, int num_atoms, /* both parameters are not used here */
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeTautString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeEquString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeIsoAtomString( INChI_IsotopicAtom   *IsotopicAtom, int nNumberOfIsotopicAtoms,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeIsoTautString( INChI_IsotopicTGroup   *IsotopicTGroup, int nNumberOfIsotopicTGroups,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeIsoHString( int num_iso_H[], char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeStereoString( AT_NUMB *at1, AT_NUMB *at2, S_CHAR *parity, int bAddDelim, int nLenCT,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeCRVString( ORIG_INFO *OrigInfo, int nLenCT, int bAddDelim,
               char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeMult( int mult, const char *szTailingDelim, char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow);
int MakeDelim( const char *szTailingDelim, char *szLinearCT, int nLen_szLinearCT, int *bOverflow);
int MakeEqStr( const char *szTailingDelim, int mult, char *szLinearCT, int nLen_szLinearCT, int *bOverflow);
int MakeHString( int bAddDelim, S_CHAR *LinearCT, int nLenCT,
                 char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow );
AT_NUMB  *GetDfsOrder4CT( AT_NUMB *LinearCT, int nLenCT, S_CHAR *nNum_H, int num_atoms, int nCtMode );

const char *EquString( int EquVal );


int str_HillFormula(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int num_components, int bUseMulipliers);
int str_HillFormula2(INCHI_SORT *pINChISort /* non-taut */, INCHI_SORT *pINChISort2 /* taut */,
                     char *pStr, int nStrLen, int tot_len,
                     int *bOverflow, int bOutType, int num_components, int bUseMulipliers);
int str_Connections(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers);
int str_H_atoms(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
               int *bOverflow, int bOutType, int ATOM_MODE, int TAUT_MODE,
               int num_components, int bUseMulipliers);
int str_Charge2(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_Sp2(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_IsoSp2(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_Sp3(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bRelRac,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_IsoSp3(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bRelRac,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_StereoAbsInv(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
               int *bOverflow, int bOutType, int num_components);
int str_IsoStereoAbsInv(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
               int *bOverflow, int bOutType, int num_components);
int str_IsoAtoms(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bAbcNumbers,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_FixedH_atoms(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers);

int str_AuxNumb(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions);
int str_AuxEqu(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_AuxTgroupEqu(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bUseMulipliers);
int str_AuxIsoTgroupEqu(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bOmitRepetitions, int bUseMulipliers);
int str_AuxInvSp3(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_AuxInvSp3Numb(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions);
int str_AuxIsoNumb(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions);
int str_AuxIsoEqu(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_AuxInvIsoSp3(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers);
int str_AuxInvIsoSp3Numb(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, char *pStr, int nStrLen, int tot_len,
              int *bOverflow, int bOutType, int TAUT_MODE, int num_components,
              int bSecondNonTautPass, int bOmitRepetitions);
int str_AuxChargeRadVal(INCHI_SORT *pINChISort, char *pStr, int nStrLen, int tot_len,
                  int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bUseMulipliers);
int bin_AuxTautTrans(INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2,
                      AT_NUMB **pTrans_n, AT_NUMB **pTrans_s, int bOutType, int num_components);
int str_AuxTautTrans(AT_NUMB *nTrans_n, AT_NUMB *nTrans_s, char *pStr, int nStrLen, int tot_len,
                     int *bOverflow, int TAUT_MODE, int num_components);
int CompareTautNonIsoPartOfINChI( const INChI *i1, const INChI *i2 );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /*__INCHIMAKE_H__*/
