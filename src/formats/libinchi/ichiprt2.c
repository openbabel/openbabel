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

#include "mode.h"

#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "ichicant.h"
#include "ichicano.h"
#include "ichicomn.h"

#include "ichicomp.h"
#include "ichimain.h"
#include "ichimake.h"

/*****************************************************************************************/
int Eql_INChI_Stereo( INChI_Stereo *s1, int eql1, INChI_Stereo *s2, int eql2, int bRelRac )
{
    int inv1=0, inv2=0, len;
    
    if ( !s1 ) {
        return 0;
    }
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
#else
    bRelRac = 0;
#endif
    
    if( EQL_SP2 == eql1 ) {
        if ( (len=s1->nNumberOfStereoBonds) > 0 && s1->b_parity && s1->nBondAtom1 && s1->nBondAtom2 ) {
            if ( !s2  ) {
                if ( EQL_EXISTS == eql2  ) {
                    /* find whether double bond stereo exists*/
                    return 1;
                }
                return 0;
            }
            if ( EQL_SP2 == eql2 &&
                 len == s2->nNumberOfStereoBonds && s2->b_parity && s2->nBondAtom1 && s2->nBondAtom2 &&
                 !memcmp( s1->nBondAtom1, s2->nBondAtom1, len * sizeof(s1->nBondAtom1[0])) &&
                 !memcmp( s1->nBondAtom2, s2->nBondAtom2, len * sizeof(s1->nBondAtom2[0])) &&
                 !memcmp( s1->b_parity,   s2->b_parity,   len * sizeof(s1->b_parity[0])) ) {
                return 1;
            }
        }
        return 0;
    } else
    if ( (eql1 == EQL_SP3 || (inv1 = (eql1 == EQL_SP3_INV))) &&
         (len=s1->nNumberOfStereoCenters) > (bRelRac? 1 : 0) ) {
        
        S_CHAR  *t_parity1, *t_parity2;
        AT_NUMB *nNumber1, *nNumber2;
        if ( inv1 ) {
            if ( s1->nCompInv2Abs ) {
                t_parity1 = s1->t_parityInv;
                nNumber1  = s1->nNumberInv;
            } else {
                return 0;
            }
        } else {
            t_parity1 = s1->t_parity;
            nNumber1  = s1->nNumber;
        }

        if ( t_parity1 && nNumber1 ) {
            if ( !s2 ) {
                if ( EQL_EXISTS == eql2 && (!inv1 || s1->nCompInv2Abs) ) {
                    /* the 1st sp3 (inverted if requested) stereo exists*/
                    return 1;
                }
                return 0;  /* both sp3 do not exist */
            }
            if( (eql2 == EQL_SP3 || (inv2 = (eql2 == EQL_SP3_INV))) &&
                len == s2->nNumberOfStereoCenters ) {
                if ( inv2 ) {
                    if ( s2->nCompInv2Abs && s1->nCompInv2Abs ) {
                        t_parity2 = s2->t_parityInv;
                        nNumber2  = s2->nNumberInv;
                    } else {
                        /* if one sp3 is inverted then another should have non-trivial inverted stereo */
                        return 0;
                    }
                } else {
                    if ( inv1 && !s2->nCompInv2Abs ) {
                        /* if one sp3 is inverted then another should have non-trivial inverted stereo */
                        return 0;
                    }
                    t_parity2 = s2->t_parity;
                    nNumber2  = s2->nNumber;
                }
                if ( t_parity2 && nNumber2 ) {
                    if ( inv1 ^ inv2 ) {
                        int i, num_inv;
                        for ( i = 0, num_inv = 0; i < len; i ++ ) {
                            if ( nNumber1[i] != nNumber2[i] )
                                break;
                            if ( ATOM_PARITY_WELL_DEF(t_parity1[i]) &&
                                 ATOM_PARITY_WELL_DEF(t_parity2[i]) ) {
                                if ( 3 == t_parity1[i] + t_parity2[i] ) {
                                    num_inv ++;
                                } else {
                                    break;
                                }
                            } else
                            if ( t_parity1[i] != t_parity2[i] ) {
                                break;
                            }
                        }
                        return (len == i && num_inv > 0); 
                    } else {
                        return !memcmp( t_parity1, t_parity2, len*sizeof(t_parity1[0])) &&
                               !memcmp( nNumber1,  nNumber2,  len*sizeof(nNumber1[0]));
                    }
                }
            }
        }
    }    
    return 0;
}
/**********************************************************************************************/
int Eql_INChI_Isotopic( INChI *i1, INChI *i2 )
{
    int eq  = i1 && i2 && !i1->bDeleted && !i2->bDeleted &&
              ( i1->nNumberOfIsotopicAtoms > 0 || i1->nNumberOfIsotopicTGroups > 0 ) &&
              i1->nNumberOfIsotopicAtoms == i2->nNumberOfIsotopicAtoms &&
              i1->nNumberOfIsotopicTGroups == i2->nNumberOfIsotopicTGroups &&
             ( !i1->nNumberOfIsotopicAtoms   ||
               i1->IsotopicAtom && i2->IsotopicAtom &&
               !memcmp( i1->IsotopicAtom,  i2->IsotopicAtom,
                        i1->nNumberOfIsotopicAtoms  * sizeof(i1->IsotopicAtom[0]) ) ) &&
             ( !i1->nNumberOfIsotopicTGroups ||
               i1->IsotopicTGroup && i2->IsotopicTGroup &&
               !memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup,
                        i1->nNumberOfIsotopicTGroups * sizeof(i1->IsotopicAtom[0]) ) );
    return eq;

}
/**********************************************************************************************/
int Eql_INChI_Aux_Equ( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 )
{
    int t1=0, t2=0, len;
    AT_NUMB *n1=NULL, *n2=NULL;
    if ( !a1 || !a2 ) {
        return 0;
    }
    t1 = (eql1 & EQL_EQU_TG);
    t2 = (eql2 & EQL_EQU_TG);
    if ( t1 && t2 ) {
        if ( (len = a1->nNumberOfTGroups) > 0 && len == a2->nNumberOfTGroups && !a1->bDeleted && !a2->bDeleted ) {
            if (eql1 & EQL_EQU_ISO) {
                if ( a1->bIsIsotopic ) {
                    n1 = a1->nConstitEquIsotopicTGroupNumbers;
                }
            } else {
                n1 = a1->nConstitEquTGroupNumbers;
            }
            if (eql2 & EQL_EQU_ISO) {
                if ( a2->bIsIsotopic ) {
                    n2 = a2->nConstitEquIsotopicTGroupNumbers;
                }
            } else {
                n2 = a2->nConstitEquTGroupNumbers;
            }
        }
    } else
    if ( !t1 && !t2 ) {
        if ( (len = a1->nNumberOfAtoms) > 0 && len == a2->nNumberOfAtoms && !a1->bDeleted && !a2->bDeleted ) {
            if (eql1 & EQL_EQU_ISO) {
                if ( a1->bIsIsotopic ) {
                    n1 = a1->nConstitEquIsotopicNumbers;
                }
            } else {
                n1 = a1->nConstitEquNumbers;
            }
            if (eql2 & EQL_EQU_ISO) {
                if ( a2->bIsIsotopic ) {
                    n2 = a2->nConstitEquIsotopicNumbers;
                }
            } else {
                n2 = a2->nConstitEquNumbers;
            }
        }
    }
    if ( n1 && n2 && !memcmp(n1, n2, len*sizeof(n1[0])) && bHasEquString( n1, len) ) {
        return 1;
    }
    return 0;
}
/**********************************************************************************************/
int Eql_INChI_Aux_Num( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 )
{
    int len;
    AT_NUMB *n1=NULL, *n2=NULL;
    if ( !a1 || !a2 ) {
        return 0;
    }
    if ( (len = a1->nNumberOfAtoms) <= 0 || len != a2->nNumberOfAtoms || a1->bDeleted || a2->bDeleted ) {
        return 0;
    }
    if ( (eql1 & EQL_NUM_ISO) && !a1->bIsIsotopic ||
         (eql2 & EQL_NUM_ISO) && !a2->bIsIsotopic ) {
        return 0;
    }

    switch ( eql1 ) {
    case EQL_NUM:
        n1 = a1->nOrigAtNosInCanonOrd;
        break;
    case EQL_NUM_ISO:
        n1 = a1->nIsotopicOrigAtNosInCanonOrd;
        break;
    case EQL_NUM_INV:
        n1 = a1->nOrigAtNosInCanonOrdInv;
        break;
    case ( EQL_NUM_INV | EQL_NUM_ISO ):
        n1 = a1->nIsotopicOrigAtNosInCanonOrdInv;
        break;
    default:
        return 0;
    }

    switch ( eql2 ) {
    case EQL_NUM:
        n2 = a2->nOrigAtNosInCanonOrd;
        break;
    case EQL_NUM_ISO:
        n2 = a2->nIsotopicOrigAtNosInCanonOrd;
        break;
    case EQL_NUM_INV:
        n2 = a2->nOrigAtNosInCanonOrdInv;
        break;
    case ( EQL_NUM_INV | EQL_NUM_ISO ):
        n2 = a2->nIsotopicOrigAtNosInCanonOrdInv;
        break;
    default:
        return 0;
    }
    
    if ( n1 && n2 && !memcmp( n1, n2, len*sizeof(n1[0])) ) {
        return 1;
    }
    return 0;
}
/**********************************************************************************************/
int bHasOrigInfo( ORIG_INFO *OrigInfo, int num_atoms )
{
    int i, bFound = 0;
    if ( OrigInfo && num_atoms > 0 ) {
        for ( i = 0; !bFound && i < num_atoms; i ++ ) {
            bFound |= (0 != OrigInfo[i].cCharge) ||
                      (0 != OrigInfo[i].cRadical) ||
                      (0 != OrigInfo[i].cUnusualValence);

        }
    }
    return bFound;
}
/**********************************************************************************************/
int EqlOrigInfo( INChI_Aux *a1, INChI_Aux *a2 )
{
    int ret =  a1 && a2 && a1->nNumberOfAtoms == a2->nNumberOfAtoms &&
               bHasOrigInfo( a1->OrigInfo, a1->nNumberOfAtoms ) && a2->OrigInfo &&
               !memcmp( a1->OrigInfo, a2->OrigInfo, a1->nNumberOfAtoms * sizeof(a1->OrigInfo[0]) );
    return ret;

}

/**********************************************************************************************/
int bHasEquString( AT_NUMB *LinearCT, int nLenCT )
{
    /*  produce output string; */
    int i, k;
    if ( !LinearCT )
        return 0;
    for ( k = 0; k < nLenCT; k ++ ) {
        /*  find the first equivalence number */
        if ( k != (int)LinearCT[k] - 1 )
            continue;
        for ( i = k; i < nLenCT; i ++ ) {
            if ( k != (int)LinearCT[i]-1 )
                continue;
            if ( k < i ) {
                return 1;
            }
        }
    }
    return 0;
}
/********************************************************************************************/
int MakeMult( int mult, const char *szTailingDelim, char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    char szValue[16];
    int  len = 0, len_delim;
    if ( mult == 1 || *bOverflow )
        return 0;
    if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
        len += MakeAbcNumber( szValue, (int)sizeof(szValue), NULL, mult );
    } else {
        len += MakeDecNumber( szValue, (int)sizeof(szValue), NULL, mult );
    }
    len_delim = strlen(szTailingDelim);
    if ( len + len_delim < (int)sizeof(szValue) ) {
        strcpy( szValue+len, szTailingDelim );
        len += len_delim;
        if ( len < nLen_szLinearCT ) {
            strcpy( szLinearCT, szValue );
            return len;
        }
    }
    *bOverflow |= 1;
    return 0;
}
/********************************************************************************************/
int MakeDelim( const char *szTailingDelim, char *szLinearCT, int nLen_szLinearCT, int *bOverflow)
{
    int len_delim;
    if ( !szTailingDelim || !*szTailingDelim || *bOverflow )
        return 0;
    len_delim = strlen(szTailingDelim);
    if ( len_delim < nLen_szLinearCT ) {
        strcpy( szLinearCT, szTailingDelim );
        return len_delim;
    }
    *bOverflow |= 1;
    return 0;
}
/********************************************************************************************/
int MakeEqStr( const char *szTailingDelim, int mult, char *szLinearCT, int nLen_szLinearCT, int *bOverflow)
{
    int  len = 0, len_delim;
    char szValue[16];
    if ( !szTailingDelim || !*szTailingDelim || *bOverflow )
        return 0;
    if ( mult != 1 ) {
        len = MakeDecNumber( szValue, (int)sizeof(szValue), NULL, mult );
    }
    len_delim = strlen(szTailingDelim);
    if ( len_delim + len < nLen_szLinearCT ) {
        if ( len > 0 ) {
            memcpy( szLinearCT, szValue, len );
        }
        strcpy( szLinearCT+len, szTailingDelim );
        return len + len_delim;
    }
    *bOverflow |= 1;
    return 0;
}
/**********************************************************************************************
 * nCtMode = 0: full
 *           1: censored CT (no orphans)
 *           2: compressed CT (Abs numbers)
 **********************************************************************************************/
int MakeCtStringNew( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  S_CHAR *nNum_H, int num_atoms,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, i, bOvfl = *bOverflow;
    char szValue[16];
    int   nValue, nDelim, num_H;
    AT_NUMB *nDfsOrderCT = NULL;
    int      bNoNum_H    = (NULL == nNum_H);
    int      nNumRingClosures;
    int      bAbcNumbers   = (0 != ( nCtMode & CT_MODE_ABC_NUMBERS ));
    int      bPredecessors = (0 != ( nCtMode & CT_MODE_PREDECESSORS ));
    int      bCountRingClosures = bAbcNumbers && bPredecessors && (nCtMode & CT_MODE_ABC_NUM_CLOSURES);
    if ( nLenCT <= 1 ) {
        return 0;  /*  no atoms or a single atom: no connection table */
    }
    /*  make array containing connection string data */
    if ( !(nDfsOrderCT = GetDfsOrder4CT( LinearCT, nLenCT, nNum_H, num_atoms, nCtMode ) ) ) {
        (*bOverflow) ++;
        return 0;
    }

    /*  add connection table string */
    if ( !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 1 ) {
            strcpy( szLinearCT, "," );
            nLen ++;
        } else {
            bOvfl = 1;
        }
    }

    if ( !bOvfl ) {
        nNumRingClosures = 0;
        for ( i = 0; nDfsOrderCT[i] && nLen < nLen_szLinearCT; i += 3 ) {
            nValue = (nDfsOrderCT[i] > MAX_ATOMS)? 0 : nDfsOrderCT[i];
            num_H  = nDfsOrderCT[i+1]? nDfsOrderCT[i+1]-16:0;
            nDelim = nDfsOrderCT[i+2];
            len = 0;
            /*  delimiter */
            if ( bPredecessors ) {
                if ( bCountRingClosures ) {
                    if ( nDelim == '-' && i > 3 && bNoNum_H ) {
                        if ( !nNumRingClosures ) {
                            int j;
                            for ( j = i; nDfsOrderCT[j] && '-' == nDfsOrderCT[j+2]; j += 3 ) {
                                nNumRingClosures ++;
                            }
                            if ( nNumRingClosures ) {
                                len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, nNumRingClosures );
                            }
                            nNumRingClosures --;
                        } else {
                            nNumRingClosures --;
                        }
                    } else {
                        nNumRingClosures = 0;
                    }
                } else
                if ( nDelim && !( bAbcNumbers && nDelim == ',' ) ) {
                    if ( nNum_H || i > 3 ) { 
                        szValue[len ++] = nDelim;
                    }
                }
            } else {
                if ( nDelim && !( bAbcNumbers && nDelim == '-' ) ) {
                    szValue[len ++] = nDelim;
                }
            }
            if ( bAbcNumbers ) {
                if ( nValue || i ) { /* the 1st value may be zero in case of presdecessor list */
                    len += MakeAbcNumber( szValue+len, (int)sizeof(szValue)-len, NULL, nValue );
                }
                if ( num_H ) {
                    len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, num_H );
                }
            } else {
                if ( nValue || i ) { /* the 1st value may be zero in case of presdecessor list */
                    len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, nValue );
                }
                if ( num_H ) {
                    szValue[len] = 'H';
                    len ++;
                    if ( num_H > 1 ) {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, num_H );
                    }
                }
            }
            if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                if ( len ) {
                    strcpy( szLinearCT+nLen, szValue );
                    nLen += len;
                }
            } else {
                bOvfl = 1;
                break;
            }
        }
    }
    *bOverflow |= bOvfl;
    if ( nDfsOrderCT )
        inchi_free( nDfsOrderCT );
    return nLen;
}
/**********************************************************************************************
 *  nCtMode = 0: full
 *            1: censored CT (no orphans)
 *            2: compressed CT (Abs numbers)
 **********************************************************************************************/
int MakeCtStringOld( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, i, bLessThanPrev, bOvfl = *bOverflow;
    AT_NUMB nMax = 0;
    char szValue[16];
    int   nValue, bNext = 0;
    /*  add connection table string */
    if ( !( nCtMode & CT_MODE_ABC_NUMBERS ) && !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 1 ) {
            strcpy( szLinearCT, "," );
            nLen ++;
        } else {
            bOvfl = 1;
        }
    }
    if ( !bOvfl ) {
        for ( i = 0; i < nLenCT && nLen < nLen_szLinearCT; i ++ ) {
            bLessThanPrev = 0;
            if ( !(nCtMode & CT_MODE_NO_ORPHANS) || ((bLessThanPrev=LinearCT[i] < nMax) ||
                  i+1 < nLenCT && LinearCT[i+1] < (nMax=LinearCT[i])) ) {
                nValue = LinearCT[i];
                if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                    len = MakeAbcNumber( szValue, (int)sizeof(szValue), (!bNext && bAddDelim)? ITEM_DELIMETER : NULL, nValue );
                } else
                if ( nCtMode & CT_MODE_NO_ORPHANS ) {  /*  censored CT */
                    /*  output '-' as a delimiter to show a bonding for decimal output of the connection table */
                    len = MakeDecNumber( szValue, (int)sizeof(szValue), bLessThanPrev? "-":ITEM_DELIMETER, nValue );
                } else {
                    len = MakeDecNumber( szValue, (int)sizeof(szValue), i? ITEM_DELIMETER:NULL, nValue );
                }
                if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                    if ( len ) {
                        strcpy( szLinearCT+nLen, szValue );
                        nLen += len;
                        bNext ++;
                    }
                } else {
                    bOvfl = 1;
                    break;
                }
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************
 *  nCtMode = 0: decimal
 *            2: compressed CT (Abs numbers)
 **********************************************************************************************/
int MakeHString( int bAddDelim, S_CHAR *LinearCT, int nLenCT,
                 char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow )
{
#define INIT_MIN_NUM_H (-4)
#define INIT_MAX_NUM_H 16
#define INIT_LEN_NUM_H (INIT_MAX_NUM_H - INIT_MIN_NUM_H + 1)

    /*  produce output string; */
    int nLen = 0, len, i, iFirst, nVal, bOvfl = *bOverflow;
    char szValue[32];
    const char *pH;
    int  bNext = 0;
    /*  add connection table string */
    if ( !( nCtMode & CT_MODE_ABC_NUMBERS ) && !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 1 ) {
            strcpy( szLinearCT, "," );
            nLen ++;
        } else {
            bOvfl = 1;
        }
    }
    if ( !bOvfl && 0 < nLenCT && LinearCT ) {
        if ( nCtMode & CT_MODE_EQL_H_TOGETHER ) {
            int  curMinH = INIT_MIN_NUM_H;
            int  curMaxH = INIT_MAX_NUM_H;
            int  curLenH = INIT_LEN_NUM_H;
            int  nInitNumH[INIT_LEN_NUM_H];
            int *nNumH = nInitNumH;
            int  numAt, curNumH;
            int      j, bOutOfRange, tot_num_no_H;
            /* count atoms H */
            do {
                bOutOfRange = 0;
                tot_num_no_H = 0; /* number of atoms that have no H */
                memset( nNumH, 0, curLenH*sizeof(nNumH[0]) );
                for ( i = 0; i < nLenCT; i ++ ) {
                    curNumH = LinearCT[i];
                    if ( curNumH < curMinH ) {
                        curMinH = curNumH;
                        bOutOfRange ++;
                    } else
                    if ( curNumH > curMaxH ) {
                        curMaxH = curNumH;
                        bOutOfRange ++;
                    } else
                    if ( !bOutOfRange ) {
                        nNumH[curNumH-curMinH] ++;
                    }
                    tot_num_no_H += !curNumH;
                }
                if ( tot_num_no_H == nLenCT ) {
                    return nLen; /* empty string */
                }
                if ( bOutOfRange ) {
                    /* for debug only */
                    if ( nNumH != nInitNumH ) {
                        *bOverflow |= 1;
                        inchi_free( nNumH );
                        return nLen;
                    }
                    /* end debug */
                    curLenH = curMaxH - curMinH + 1;
                    nNumH = (int*) inchi_malloc( curLenH * sizeof(nNumH[0]) );
                    if ( !nNumH ) {
                        *bOverflow |= 1;
                        return nLen;
                    }
                }
            } while ( bOutOfRange ); /* the loop may be executed 1 or 2 times only */

            for ( curNumH = curMinH;  curNumH <= curMaxH; curNumH ++ ) {
                numAt = nNumH[curNumH-curMinH]; /* number of atoms that have curNumH atoms H */
                if ( !numAt || !curNumH ) {
                    continue; /* no atom has this number of H or number of H = 0 */
                }
                j = 0;
                while ( j < nLenCT && numAt ) {
                    if ( curNumH == LinearCT[j] ) {
                        iFirst = ++j;
                        numAt --;
                        for ( ; j < nLenCT && curNumH == LinearCT[j] && numAt; j ++ ) {
                            numAt --;
                        }
                        if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                            len = MakeAbcNumber( szValue, (int)sizeof(szValue), NULL, iFirst );
                        } else {
                            len = MakeDecNumber( szValue, (int)sizeof(szValue), bNext?ITEM_DELIMETER:NULL, iFirst );
                            bNext ++; /* add a delimiter (comma) before all except the first */
                        }
                        if ( iFirst < j ) {
                            /* output last canonical number */
                            if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                                len += MakeAbcNumber( szValue+len, (int)sizeof(szValue), NULL, j );
                            } else {
                                len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, "-", j );
                            }
                        }
                        if ( !numAt || ( nCtMode & CT_MODE_ABC_NUMBERS ) ) {
                            /* add number of H */
                            /* output number of H */
                            nVal = curNumH;
                            if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                                len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, nVal );
                            } else {
                                pH = nVal > 0? "H":"h";
                                nVal = abs(nVal);
                                if ( nVal > 1 ) {
                                    len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, pH, nVal );
                                } else {
                                    strcpy( szValue+len, pH );
                                    len ++;
                                }
                            }
                        }
                        /* add to the output */
                        if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                            if ( len ) {
                                strcpy( szLinearCT+nLen, szValue );
                                nLen += len;
                                bNext ++;
                            }
                        } else {
                            bOvfl = 1;
                            break;
                        }
                    } else {
                        j ++;
                    }
                }
            }
            if ( nNumH != nInitNumH ) {
                inchi_free( nNumH );
            }
        } else {
            iFirst     = 0;
            for ( i = iFirst+1; i <= nLenCT && nLen < nLen_szLinearCT; i ++ ) {
                if ( i < nLenCT && LinearCT[i] == LinearCT[iFirst] ) {
                    continue;
                }
                /* output identical values located at i = iFirst..i-1 */
                if ( LinearCT[iFirst] ) { /* output only non-zero values */
                    /* first canonical number */
                    nVal = LinearCT[iFirst];
                    iFirst ++;
                    if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                        len = MakeAbcNumber( szValue, (int)sizeof(szValue), NULL, iFirst );
                    } else {
                        len = MakeDecNumber( szValue, (int)sizeof(szValue), bNext?ITEM_DELIMETER:NULL, iFirst );
                    }
                    if ( iFirst < i ) {
                        /* output last canonical number */
                        if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                            len += MakeAbcNumber( szValue+len, (int)sizeof(szValue), NULL, i );
                        } else {
                            len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, "-", i );
                        }
                    }
                    /* output number of H */
                    if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, nVal );
                    } else {
                        pH = nVal > 0? "H":"h";
                        nVal = abs(nVal);
                        if ( nVal > 1 ) {
                            len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, pH, nVal );
                        } else {
                            strcpy( szValue+len, pH );
                            len ++;
                        }
                    }
                    if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                        if ( len ) {
                            strcpy( szLinearCT+nLen, szValue );
                            nLen += len;
                            bNext ++;
                        }
                    } else {
                        bOvfl = 1;
                        break;
                    }
                }
                iFirst = i;
            }
        }
    }

    *bOverflow |= bOvfl;
    return nLen;

#undef INIT_MIN_NUM_H
#undef INIT_MAX_NUM_H
#undef INIT_LEN_NUM_H
}
/**********************************************************************************************
 *  nCtMode = 0: full
 *            1: censored CT (no orphans, that CT should have only atoms with neighbors)
 *            2: compressed CT (Abc numbers)
 **********************************************************************************************/
int MakeCtString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  S_CHAR *nNum_H, int num_atoms, /* both parameters are not used here */
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    
    if ( !nNum_H || !(nCtMode & CT_MODE_NO_ORPHANS) ) {
        return MakeCtStringOld( LinearCT, nLenCT, bAddDelim,
                                szLinearCT, nLen_szLinearCT, nCtMode, bOverflow);
    } else {
        return MakeCtStringNew( LinearCT, nLenCT, bAddDelim,
                             nNum_H, num_atoms,
                             szLinearCT, nLen_szLinearCT, nCtMode, bOverflow);
    }
}    



/**********************************************************************************************
 *  nCtMode = 0: full: decimal-only, with parentheses around t-groups
 *            2: compressed CT: do not add comma before the output string if bAddDelim != 0
 *                            do not add parentheses around t-groups
 *                            atom canon numbers an Abc
 * LinearCT format:
 *    N      = number of tautomeric groups
 *    n      = number of endpoints + 1 in a tautomeric group #1
 *   next INCHI_T_NUM_MOVABLE lines (any after the first non-zero):
 *    h      = number of hydrogen atoms in the tautomeric group
 *    m      = number of negative charges
 *    ...    (the rest of the INCHI_T_NUM_MOVABLE has not been established, ignore them)
 *    c(1)   = canonical number of the first atom in the t-group
 *    ...
 *    c(n-1) = canonical number of the last atom in the t-group
 *
 **********************************************************************************************/

int MakeTautString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, i, bOvfl = *bOverflow;
    char szValue[16];
    const char *p;
    int   nValue, nGroupLen, iGroupOutputCount, bCompressed;
    /*  make tautomer string */
    if ( !nLenCT || !LinearCT || !*LinearCT ) {
        return nLen;
    }
    bCompressed = ( nCtMode & CT_MODE_ABC_NUMBERS );
    if ( !bCompressed && !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 1+LEN_EXTRA_SPACE ) {
            strcpy( szLinearCT, COMMA_EXTRA_SPACE);
            nLen += 1+LEN_EXTRA_SPACE;
        } else {
            bOvfl = 1;
        }
    }
    LinearCT ++; /*  bypass number of tautomeric groups */
    nLenCT --;
    
    if ( !bOvfl ) {
        for ( i = nGroupLen = iGroupOutputCount = 0; i < nLenCT && nLen < nLen_szLinearCT; i ++ ) {
            nValue = (int)LinearCT[i];
            if ( nGroupLen == iGroupOutputCount ) {
                nGroupLen = nValue;
                iGroupOutputCount = 0;
                /* group delimiter (uncompressed) */
                if ( !bCompressed ) {
                    if ( !i ) {
                        strcpy( szValue, "(" );
                        len = 1;
                    } else {
                        strcpy( szValue, ")(" );
                        len = 2;
                    }
                } else {
                    len = 0;
                }
            } else
            if ( bCompressed && iGroupOutputCount >= INCHI_T_NUM_MOVABLE ) {
                /*  compressed canon number in Abc */
                len = MakeAbcNumber( szValue, (int)sizeof(szValue), NULL, nValue );
                iGroupOutputCount ++;
            } else {
                /*  always output number of hydrogen atoms as a decimal */
                /*  output leading space if: */
                /*  (a) this is the first output value in compressed mode (i==1 && bCompressed) */
                /*  (b) this is not the first output value in non-compressed mode ( iGroupOutputCount && !bCompressed) */
                if ( bCompressed ) {
                    p   = NULL;
                    len = 0;
                    switch( iGroupOutputCount ) {
                    case 0:
                        len = MakeDecNumber( szValue, (int)sizeof(szValue), (i == 1)? ITEM_DELIMETER:NULL, nValue );
                        break;
                    case 1:
                        p = "-";
                        break;
                    case 2:
                        p = "+";
                        break;
                    }
                    if ( p ) {
                        switch( nValue ) {
                        case 0:
                            len = 0;
                            break;
                        case 1:
                            strcpy(szValue, p);
                            len = strlen(szValue);
                            break;
                        default:
                            len = MakeDecNumber( szValue, (int)sizeof(szValue), p, nValue );
                            break;
                        }
                    }
                } else {
                    if ( iGroupOutputCount >= INCHI_T_NUM_MOVABLE ) {
                        /*  canonical number of the atom in the tautomeric group */
                        len = MakeDecNumber( szValue, (int)sizeof(szValue), ITEM_DELIMETER, nValue );
                    } else {
                        p   = NULL;
                        len = 0;
                        if ( nValue ) {
                            switch( iGroupOutputCount ) {
                            case 0:
                                p = "H";
                                break;
                            case 1:
                                p = "-";
                                break;
                            case 2:
                                p = "+";
                                break;
                            }
                            if ( p ) {
                                /*  number of hydrogens */
                                if ( nValue == 1 ) {
                                    strcpy(szValue, p);
                                    len = strlen(szValue);
                                } else {
                                    len = MakeDecNumber( szValue, (int)sizeof(szValue), p, nValue );
                                }
                            }
                        }
                    }
                }
                iGroupOutputCount ++;
            }
            if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                if ( len ) {
                    strcpy( szLinearCT+nLen, szValue );
                    nLen += len;
                }
            } else {
                bOvfl = 1;
                break;
            }
        }
        if ( !bOvfl && !bCompressed && i ) {
            if ( nLen + 1 < nLen_szLinearCT ) {
                strcpy( szLinearCT+nLen, ")" );
                nLen ++;
            } else {
                bOvfl = 1;
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************
 *  nCtMode = 0: full
 *            2: compressed CT
 *  22+3s3:  22=canon. number; +3=charge; s=singlet (d=doublet, t=triplet, s is omitted if valence=0), 3 = valence
 *  22+3.3, (charge, valence) 22.3 (valence) 22t3 (triplet, valence) 
 *  Ab+3t4:  Ab=canon. number; +3=charge or "." t=triplet (or s, d), 4=valence
 **********************************************************************************************/
int MakeCRVString( ORIG_INFO *OrigInfo, int nLenCT, int bAddDelim,
               char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, k, bAbcNumbers;
    int bOvfl = *bOverflow;
    char szValue[32];
    int  bNext=0;
    bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    /*  add connection table string */
    if ( !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 2 ) {
            strcpy( szLinearCT, ", " );
            nLen += 2;
        } else {
            bOvfl = 1;
        }
    }
    for ( k = 0; !bOvfl && k < nLenCT && nLen < nLen_szLinearCT; k ++ ) {
        /*  find the next non-empty entry */
        if ( OrigInfo[k].cCharge || OrigInfo[k].cRadical || OrigInfo[k].cUnusualValence ) {
            if ( bAbcNumbers ) {
                /*
                    3 items: Ad+3d4 (canon. numb=Ad, charge=+3, doublet, valence = 4
                    2 items: Ad.d4  Ad+3.4  Ad+3d
                    1 item:  Ad+3   Ad.d     Ad4

                    dot output before radical: no charge, radical is present
                    dot before valence:        charge is present, no radical, valence is present
                 */
                len = MakeAbcNumber( szValue, (int)sizeof(szValue), NULL, k+1 );

                /* charge */
                if ( OrigInfo[k].cCharge ) {
                    if ( OrigInfo[k].cCharge > 0 ) {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, "+", OrigInfo[k].cCharge );
                    } else {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, OrigInfo[k].cCharge );
                    }
                }
                /* radical */
                if ( OrigInfo[k].cRadical ) {
                    if ( !OrigInfo[k].cCharge ) {
                        szValue[len ++] = '.';
                    }
                    switch( OrigInfo[k].cRadical ) {
                    case 1:
                        szValue[len ++] = 'd';
                        break;
                    case 2:
                        szValue[len ++] = 't';
                        break;
                    default:
                        szValue[len ++] = 'u';
                        break;
                    }
                }
                /* valence */
                if ( OrigInfo[k].cUnusualValence ) {
                    if ( OrigInfo[k].cCharge && !OrigInfo[k].cRadical ) {
                        szValue[len ++] = '.';
                    }
                    len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, OrigInfo[k].cUnusualValence );
                }
            } else {
                /*
                    3 items: 22+3d4 (canon. numb=22, charge=+3, doublet, valence = 4
                    2 items: 22d4  22+3.4  22+3d
                    1 item:  22+3  22d     22.4

                    dot output before valence:
                                (a) charge,    no radical, valence
                                (b) no charge, no radical, valence
                    that is, whenever valence is present and no radical
                 */
                len = MakeDecNumber( szValue, (int)sizeof(szValue), bNext? ITEM_DELIMETER:NULL, k+1 );
                /* charge */
                if ( OrigInfo[k].cCharge ) {
                    if ( OrigInfo[k].cCharge > 0 ) {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, "+", OrigInfo[k].cCharge );
                    } else {
                        len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, OrigInfo[k].cCharge );
                    }
                }
                /* radical */
                if ( OrigInfo[k].cRadical ) {
                    switch( OrigInfo[k].cRadical ) {
                    case 1:
                        szValue[len ++] = 'd';
                        break;
                    case 2:
                        szValue[len ++] = 't';
                        break;
                    default:
                        szValue[len ++] = 'u';
                        break;
                    }
                }
                /* valence */
                if ( OrigInfo[k].cUnusualValence ) {
                    if ( !OrigInfo[k].cRadical ) {
                        szValue[len ++] = '.';
                    }
                    len += MakeDecNumber( szValue+len, (int)sizeof(szValue)-len, NULL, OrigInfo[k].cUnusualValence );
                }
            }
        } else {
            len = 0;
        }
        if ( len && nLen+len < nLen_szLinearCT ) {
            strcpy( szLinearCT+nLen, szValue );
            nLen += len;
            bNext ++;
        } else
        if ( len ) {
            bOvfl = 1;
            break;
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}

/**********************************************************************************************
 *  nCtMode = 0: full
 *            2: compressed CT
 **********************************************************************************************/
int MakeEquString( AT_NUMB *LinearCT, int nLenCT, int bAddDelim,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, i, k, bAbcNumbers;
    int bOvfl = *bOverflow;
    char szValue[16];
    int  bNext=0;
    bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    /*  add connection table string */
    if ( !bOvfl && bAddDelim ) {
        if ( nLen_szLinearCT > 2 ) {
            strcpy( szLinearCT, ", " );
            nLen += 2;
        } else {
            bOvfl = 1;
        }
    }
    for ( k = 0; !bOvfl && k < nLenCT && nLen < nLen_szLinearCT; k ++ ) {
        /*  find the first equivalence number */
        if ( k != (int)LinearCT[k] - 1 )
            continue;
        for ( i = k; i < nLenCT && nLen < nLen_szLinearCT; i ++ ) {
            if ( k != (int)LinearCT[i]-1 )
                continue;
            /*  equivalence number: a minimal canon_number out of a group of equivalent atoms */
            /*  is at canon_number-1 position of each equivalent atom.  */
            if ( bAbcNumbers ) {
                len = MakeAbcNumber( szValue, (int)sizeof(szValue), (i==k && bNext)? ITEM_DELIMETER : NULL, i+1 );
            } else {
                len = MakeDecNumber( szValue, (int)sizeof(szValue), (i==k)? "(":ITEM_DELIMETER, i+1 );
            }
            if ( 0 <= len && nLen+len < nLen_szLinearCT ) {
                strcpy( szLinearCT+nLen, szValue );
                nLen += len;
                bNext ++;
            } else
            if ( 0 > len ) {
                bOvfl = 1;
                break;
            }
        }
        if ( !bOvfl && !bAbcNumbers ) {
            if ( nLen + 2 < nLen_szLinearCT ) {
                strcpy( szLinearCT+nLen, ")" );
                nLen ++;
            } else {
                bOvfl = 1;
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************
 *  nCtMode = 0: full
 *            2: compressed CT
 **********************************************************************************************/
int MakeIsoAtomString( INChI_IsotopicAtom   *IsotopicAtom, int nNumberOfIsotopicAtoms,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, tot_len, ret, i, j, bOvfl = *bOverflow;
    char szValue[64];
    char *p;
    int   nValue;
    int   bAbcNumbers = (nCtMode & CT_MODE_ABC_NUMBERS );
    static const char letter[]  = "itdh";
    static const char *h[]      = {"T", "D", "H"};
    static const char *sign[]   = {"-", "+"};

    if ( !bOvfl ) {
        for ( i = 0; i < nNumberOfIsotopicAtoms && nLen < nLen_szLinearCT; i ++ ) {
            p = szValue;
            tot_len = 0;
            for ( j = 0; j < 5; j ++ ) {
                len = 0;
                switch( j ) {
                case 0:
                    nValue = (int)IsotopicAtom[i].nAtomNumber;
                    break;
                case 1:
                    nValue = (int)IsotopicAtom[i].nIsoDifference;
                    break;
                case 2:
                    nValue = (int)IsotopicAtom[i].nNum_T;
                    break;
                case 3:
                    nValue = (int)IsotopicAtom[i].nNum_D;
                    break;
                case 4:
                    nValue = (int)IsotopicAtom[i].nNum_H;
                    break;
                }
                if ( !j ) {
                    /*  atom canonical number */
                    len = (bAbcNumbers? MakeAbcNumber:MakeDecNumber)
                               ( p, (int)sizeof(szValue)-tot_len,
                                 bAbcNumbers?NULL:(i?ITEM_DELIMETER:EXTRA_SPACE), nValue
                               );
                } else
                if ( bAbcNumbers ) { /*  Abc output */
                    switch ( j ) {
                    case 1: /* nIsoDifference */
                        len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, NULL, nValue );
                        break;
                    case 2: /* nNum_T */
                    case 3: /* nNum_D */
                    case 4: /* nNum_H */
                        if ( nValue ) {
                            if ( (int)sizeof(szValue) - tot_len > 1 ) {
                                p[len++]=letter[j-1];
                                if ( 1 == nValue ) {
                                    p[len]  = '\0';
                                } else {
                                    ret = MakeDecNumber( p+len, (int)sizeof(szValue)-tot_len-len, NULL, nValue );
                                    len = (ret >= 0)? len+ret : ret;
                                }
                            } else {
                                len = -1; /* overflow */
                            }
                        }
                    }
                } else
                if ( nValue ) {
                    if ( j == 1 ) { /*  Decimal output */
                        /*  signed isotopic mass difference */
                        int subtract = (nValue > 0);
                        /*  (n = mass difference) > 0 corresponds to nValue = n+1 */
                        /*  subtract 1 from it so that mass difference for 35Cl or 12C is zero */
                        len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, sign[nValue>=0], abs(nValue-subtract) );
                    } else {
                        /*  hydrogen isotope */
                        if ( nValue != 1 ) {
                            len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, h[j-2], nValue );
                        } else
                        if ( (int)sizeof(szValue)-tot_len > 1 ) {
                            strcpy( p, h[j-2] );
                            len = 1;
                        } else {
                            len = -1; /*  overflow */
                        }
                    }
                } else {
                    continue; /*  do not write zeroes */
                }
                if ( len < 0 ) {
                    bOvfl = 1;
                    break;
                }
                tot_len += len;
                p += len;
            }
            if ( nLen+tot_len < nLen_szLinearCT ) {
                memcpy( szLinearCT+nLen, szValue, tot_len+1 );
                nLen += tot_len;
            } else {
                bOvfl = 1;
                break;
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************/
int MakeIsoTautString( INChI_IsotopicTGroup   *IsotopicTGroup, int nNumberOfIsotopicTGroups,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, tot_len, i, j, bOvfl = *bOverflow;
    AT_NUMB nMax;
    char szValue[32];
    char *p;
    int   nValue;
    int   bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    static const char letter[] = "tdh";
    static const char *h[]     = {"T", "D", "H"};
    /*  add connection table string */
    nMax = 0;
    if ( !bOvfl ) {
        for ( i = 0; i < nNumberOfIsotopicTGroups && nLen < nLen_szLinearCT; i ++ ) {
            p = szValue;
            tot_len = 0;
            for ( j = 0; j < 4; j ++ ) {
                switch( j ) {
                case 0:
                    nValue = (int)IsotopicTGroup[i].nTGroupNumber;
                    break;
                case 1:
                    nValue = (int)IsotopicTGroup[i].nNum_T;
                    break;
                case 2:
                    nValue = (int)IsotopicTGroup[i].nNum_D;
                    break;
                case 3:
                    nValue = (int)IsotopicTGroup[i].nNum_H;
                    break;
                }
                if ( !j ) {
                    /*  atom canonical number */
                    len = (bAbcNumbers?MakeAbcNumber:MakeDecNumber)
                            ( p, (int)sizeof(szValue)-tot_len,
                              bAbcNumbers?NULL:(i?ITEM_DELIMETER:EXTRA_SPACE),
                              nValue
                            );
                } else
                if ( nValue ) {
                    if ( bAbcNumbers ) {
                        len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, NULL, nValue );
                        if ( len > 0 ) { /*  make sure overflow has not happened */
                            if ( (int)sizeof(szValue)-tot_len-len > 1 ) {
                                p[len++]=letter[j-1];
                                p[len]  = '\0';
                            } else {
                                len = -1; /*  overflow */
                            }
                        }
                    } else {
                        /*  hydrogen isotope */
                        if ( nValue != 1 ) {
                            len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, h[j-1], nValue );
                        } else
                        if ( (int)sizeof(szValue)-tot_len > 1 ) {
                            strcpy( p, h[j-1] );
                            len = 1;
                        } else {
                            len = -1; /*  overflow */
                        }
                    }
                } else {
                    continue; /*  do not write zeroes */
                }
                if ( len < 0 ) {
                    bOvfl = 1;
                    break;
                }
                p += len;
                tot_len += len;
            }
            if ( nLen+tot_len < nLen_szLinearCT ) {
                memcpy( szLinearCT+nLen, szValue, tot_len+1 );
                nLen += tot_len;
            } else {
                bOvfl = 1;
                break;
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************/
int MakeIsoHString( int num_iso_H[], char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, tot_len, j, bOvfl = *bOverflow;
    AT_NUMB nMax;
    char szValue[32];
    char *p;
    int   nValue;
    int   bAbcNumbers = ( nCtMode & CT_MODE_ABC_NUMBERS );
    static const char letter[] = "tdh";
    static const char *h[]     = {"T", "D", "H"};
    /*  add connection table string */
    nMax = 0;
    if ( !bOvfl ) {
        p = szValue;
        tot_len = 0;
        for ( j = 1; j < 4; j ++ ) {
            nValue = num_iso_H[NUM_H_ISOTOPES-j];/* j: 1=>T, 2=>D, 3=>1H */
            if ( nValue ) {
                if ( bAbcNumbers ) {
                    len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, NULL, nValue );
                    if ( len > 0 ) { /*  make sure overflow has not happened */
                        if ( (int)sizeof(szValue)-tot_len-len > 1 ) {
                            p[len++]=letter[j-1];
                            p[len]  = '\0';
                        } else {
                            len = -1; /*  overflow */
                        }
                    }
                } else {
                    /*  hydrogen isotope */
                    if ( nValue != 1 ) {
                        len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, h[j-1], nValue );
                    } else
                    if ( (int)sizeof(szValue)-tot_len > 1 ) {
                        strcpy( p, h[j-1] );
                        len = 1;
                    } else {
                        len = -1; /*  overflow */
                    }
                }
            } else {
                continue; /*  do not write zeroes */
            }
            if ( len < 0 ) {
                bOvfl = 1;
                break;
            }
            p += len;
            tot_len += len;
        }
        if ( nLen+tot_len < nLen_szLinearCT ) {
            memcpy( szLinearCT+nLen, szValue, tot_len+1 );
            nLen += tot_len;
        } else {
            bOvfl = 1;
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
/**********************************************************************************************/
int MakeStereoString( AT_NUMB *at1, AT_NUMB *at2, S_CHAR *parity, int bAddDelim, int nLenCT,
                  char *szLinearCT, int nLen_szLinearCT, int nCtMode, int *bOverflow)
{
    /*  produce output string; */
    int nLen = 0, len, tot_len, i, j, bOvfl = *bOverflow;
    char szValue[32];
    char *p;
    int   nValue;
    static const char parity_char[] = "!-+u?";
    bAddDelim = 0;
    if ( !bOvfl ) {
        for ( i = 0; i < nLenCT && nLen < nLen_szLinearCT; i ++ ) {
            p = szValue;
            tot_len = 0;
            for ( j = 0; j < 3; j ++ ) {
                if ( j == 0 && at1 )
                    nValue = (int)at1[i];
                else
                if ( j == 1 && at2 )
                    nValue = (int)at2[i];
                else
                if ( j == 2 && parity )
                    nValue = (int)parity[i];
                else
                    continue;
                if ( nCtMode & CT_MODE_ABC_NUMBERS ) {
                    len = (j==2? MakeDecNumber : MakeAbcNumber)( p, (int)sizeof(szValue)-tot_len, NULL, nValue );
                } else {
                    if ( j < 2 ) {
                        len = MakeDecNumber( p, (int)sizeof(szValue)-tot_len, tot_len?"-":(i||bAddDelim)?ITEM_DELIMETER:NULL, nValue );
                    } else
                    if ( tot_len + 1 < (int)sizeof(szValue) ) {
                        *p ++ = (0<=nValue && nValue<=4)? parity_char[nValue]:parity_char[0];
                        *p    = '\0';
                        len   = 1;
                    } else {
                        len = -1; /*  Overflow */
                    }
                }
                if ( len < 0 ) {
                    bOvfl = 1;
                    break;
                }
                p += len;
                tot_len += len;
            }
            if ( nLen+tot_len < nLen_szLinearCT ) {
                memcpy( szLinearCT+nLen, szValue, tot_len+1 );
                nLen += tot_len;
            } else {
                bOvfl = 1;
                break;
            }
        }
    }
    *bOverflow |= bOvfl;
    return nLen;
}
#ifdef ALPHA_BASE

#if ( ALPHA_BASE != 27 )
#error ALPHA_BASE definitions mismatch
#endif

#else

#define ALPHA_BASE     27

#endif

#define ALPHA_MINUS    '-'
#define ALPHA_ZERO_VAL '.'
#define ALPHA_ONE      'a'
#define ALPHA_ZERO     '@'
/**********************************************************************************************/
/*  Produce an "Alphabetic" number, base 27 (27 digits: 0, a, b, ..., z) */
/*  The leading "digit" uppercase, the rest -- lowercase */
/*  szString length nStringLen includes 1 byte for zero termination */
/*  Return Value: length without zero termination; -1 means not enough room */
/*  Note: ASCII-encoding specific implementation */
int MakeAbcNumber( char *szString, int nStringLen, const char *szLeadingDelim, int nValue )
{
    char *p = szString;
    char *q;
    int  nChar;

    if ( nStringLen < 2 )
        return -1;
    while ( szLeadingDelim && *szLeadingDelim && --nStringLen ) {
        *p ++ = *szLeadingDelim ++;
    }
    if ( nStringLen < 2 )
        return -1;
    if ( !nValue ) {
        *p++ = ALPHA_ZERO_VAL;  /*  zero value (cannot use 0) */
        *p   = '\0';
        return 1;
    }
    if ( nValue < 0 ) {
        *p++ = ALPHA_MINUS;
        nStringLen --;
        nValue = -nValue;
    }
    for ( q = p; nValue && --nStringLen; nValue /= ALPHA_BASE ) {
        if ( nChar = nValue % ALPHA_BASE ) {
            nChar = ALPHA_ONE + nChar - 1;
        } else {
            nChar = ALPHA_ZERO;
        }
        *q++ = nChar;
    }
    if ( nStringLen <= 0 )
        return -1;
    *q = '\0';
    mystrrev( p );
    p[0] = toupper(p[0]);
    return (q - szString);
}
#if ( READ_INCHI_STRING == 1 )
/*****************************************************/
static long abctol( const char *szString, char **q ); /* keep compiler happy */

long abctol( const char *szString, char **q )
{
#define __MYTOLOWER(c) ( ((c) >= 'A') && ((c) <= 'Z') ? ((c) - 'A' + 'a') : (c) )

    long        val  = 0;
    long        sign = 1;
    const char *p = szString;
    if ( *p == ALPHA_MINUS ) {
        p ++;
        sign = -1;
    }
    if ( *p == ALPHA_ZERO ) {
        p ++;
        goto exit_function;
    }
    if ( !isupper(UCINT *p) ) {
        p = szString;
        goto exit_function; /* not an abc-number */
    }
    val = __MYTOLOWER(*p) - ALPHA_ONE + 1;
    p ++;
    while ( *p ) {
        if ( islower( UCINT *p ) ) {
            val *= ALPHA_BASE;
            val += *p - ALPHA_ONE + 1;
        } else
        if ( *p == ALPHA_ZERO ) {
            val *= ALPHA_BASE;
        } else {
            break;
        }
        p ++;
    }
exit_function:
    if ( q ) {
        *q = (char *)p;  /* cast deliberately discards const qualifier */
    }
    return val;
#undef __MYTOLOWER
}
/********************************************************/
long inchi_strtol( const char *str, const char **p, int base)
{
    if ( base == ALPHA_BASE ) {
        return abctol( str, (char **)p ); /* cast deliberately discards const qualifier */
    } else {
        return strtol( str, (char **)p, base ); /* cast deliberately discards const qualifier */
    }
}
#endif
#undef ALPHA_BASE
#undef ALPHA_MINUS
#undef ALPHA_ZERO_VAL
#undef ALPHA_ONE
#undef ALPHA_ZERO

/********************************************************/
double inchi_strtod( const char *str, const char **p )
{
        return strtod( str, (char **)p );
}

/**********************************************************************************************/
/*  Produce a decimal number */
/*  szString length nStringLen includes 1 byte for zero termination */
/*  Return Value: length without zero termination; -1 means not enough room */
int MakeDecNumber( char *szString, int nStringLen, const char *szLeadingDelim, int nValue )
{
#define DECIMAL_BASE     10
#define DECIMAL_MINUS    '-'
#define DECIMAL_ZERO_VAL '0'
#define DECIMAL_ONE      '1'
#define DECIMAL_ZERO     '0'
    char *p = szString;
    char *q;
    int  nChar;

    if ( nStringLen < 2 )
        return -1;
    while ( szLeadingDelim && *szLeadingDelim && --nStringLen ) {
        *p ++ = *szLeadingDelim ++;
    }
    if ( nStringLen < 2 )
        return -1;
    if ( !nValue ) {
        *p++ = DECIMAL_ZERO_VAL;  /*  zero value (cannot use 0) */
        *p   = '\0';
        return p-szString;
    }
    if ( nValue < 0 ) {
        *p++ = DECIMAL_MINUS;
        nStringLen --;
        nValue = -nValue;
    }
    for ( q = p; nValue && --nStringLen; nValue /= DECIMAL_BASE ) {
        if ( nChar = nValue % DECIMAL_BASE ) {
            nChar = DECIMAL_ONE + nChar - 1;
        } else {
            nChar = DECIMAL_ZERO;
        }
        *q++ = nChar;
    }
    if ( nStringLen <= 0 )
        return -1;
    *q = '\0';
    mystrrev( p );
    return (q - szString);
#undef DECIMAL_BASE
#undef DECIMAL_MINUS
#undef DECIMAL_ZERO_VAL
#undef DECIMAL_ONE
#undef DECIMAL_ZERO
}
