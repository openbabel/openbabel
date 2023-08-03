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
#include <math.h>


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


int GetHillFormulaCounts( U_CHAR *nAtom, S_CHAR *nNum_H, int num_atoms,
                          AT_NUMB *nTautomer, int lenTautomer,
                          int *pnum_C, int *pnum_H, int *pnLen, int *pnNumNonHAtoms );
int MakeHillFormula( U_CHAR *nAtom, int num_atoms,
                  char *szLinearCT, int nLen_szLinearCT, int num_C, int num_H, int *bOverflow );

#if ( FIX_DALKE_BUGS == 1 )
#else
char *AllocateAndFillHillFormula( INChI *pINChI );
#endif

int AddElementAndCount( const char *szElement, int mult, char *szLinearCT, int nLenLinearCT, int *bOverflow );

int Copy2StereoBondOrAllene( INChI_Stereo *Stereo, int *nNumberOfStereoCenters, int *nNumberOfStereoBonds,
                            AT_STEREO_DBLE *LinearCTStereoDble,
                            AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic );

int CopyLinearCTStereoToINChIStereo( INChI_Stereo *Stereo,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble
           , AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic
           , AT_STEREO_CARB *LinearCTStereoCarbInv
           , AT_STEREO_DBLE *LinearCTStereoDbleInv
           , AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv );
int GetHillFormulaIndexLength( int count );

int MarkAmbiguousStereo( sp_ATOM *at, inp_ATOM *norm_at, int bIsotopic, AT_NUMB *pCanonOrd,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble );

INCHI_MODE UnmarkAllUndefinedUnknownStereo( INChI_Stereo *Stereo, INCHI_MODE nUserMode );

int CleanCoord( MOL_COORD szCoord, int delim );

/**********************************************************************************************/
int MakeHillFormulaString( char *szHillFormula, char *szLinearCT, int nLen_szLinearCT, int *bOverflow)
{
    int nLen;
    if ( szHillFormula && !*bOverflow ) {
        if ( nLen_szLinearCT > ( nLen = strlen(szHillFormula) ) ) {
            memcpy( szLinearCT, szHillFormula, nLen+1 );
            return nLen;
        }
        *bOverflow |= 1;
        return nLen_szLinearCT+1;
    }
    return 0;
}
/**********************************************************************************************
 *     MS Windows dependent: sprintf() is supposed to return the length of the output string
 *     Carbon atoms are always first
 *     Bridging hydrogen atoms are always last
 **********************************************************************************************/
int GetHillFormulaIndexLength( int count )
{
    char szCount[16];
    if ( count > 1 ) {
        return sprintf( szCount, "%d", count );
    }
    return 0;
}
/**********************************************************************************************/
int GetHillFormulaCounts( U_CHAR *nAtom, S_CHAR *nNum_H, int num_atoms,
                          AT_NUMB *nTautomer, int lenTautomer,
                          int *pnum_C, int *pnum_H, int *pnLen, int *pnNumNonHAtoms )
{
    char szElement[4];
    U_CHAR nPrevAtom = (U_CHAR)-2;
    int  bCarbon, bHydrogen, nElemLen, nFormLen, nNumNonHAtoms;
    int  mult, i, num_H, num_C;

    num_H     = 0;
    num_C     = 0;
    bCarbon   = 0;
    bHydrogen = 0;
    nElemLen  = 0;
    nFormLen  = 0;
    mult      = 0;
    nNumNonHAtoms = num_atoms;
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( nPrevAtom != nAtom[i] ) {
            if ( mult ) {
                if ( bHydrogen ) {
                    num_H += mult;
                }else
                if ( bCarbon ) {
                    num_C += mult;
                } else {
                    nFormLen += nElemLen;
                    nFormLen += GetHillFormulaIndexLength( mult );
                }
            }

            if ( GetElementFormulaFromAtNum((int)nAtom[i], szElement ) ) {
                return -1; /*  wrong element */
            }
            mult = 1;

            nElemLen = strlen(szElement);
            nPrevAtom = nAtom[i];
            bCarbon   = !strcmp( szElement, "C" );
            bHydrogen = !strcmp( szElement, "H" );
            if ( bHydrogen ) {
                nNumNonHAtoms = i;
            }
        } else {
            mult ++;
        }
        
        num_H += nNum_H[i];
    }
    /* NumGroups; ((NumAt+1, NumH, At1..AtNumAt),...) */
    if ( nTautomer && lenTautomer > 0 ) {
        int num_groups = nTautomer[0];
        for ( i = 1; i < lenTautomer && num_groups > 0; i += nTautomer[i]+1, num_groups -- ) {
            num_H += nTautomer[i+1];
        }
    }

    if ( mult ) {
        if ( bHydrogen ) {
            num_H += mult;
        } else
        if ( bCarbon ) {
            num_C += mult;
        } else {
            nFormLen += nElemLen;
            nFormLen += GetHillFormulaIndexLength( mult );
        }
    }

    if ( num_C ) {
        nFormLen += strlen( "C" );
        nFormLen += GetHillFormulaIndexLength( num_C );
    }

    if ( num_H ) {
        nFormLen += strlen( "H" );
        nFormLen += GetHillFormulaIndexLength( num_H );
    }
    *pnum_C = num_C;
    *pnum_H = num_H;
    *pnLen  = nFormLen;
    *pnNumNonHAtoms = nNumNonHAtoms;

    return 0;
}
/**********************************************************************************************/
int AddElementAndCount( const char *szElement, int mult, char *szLinearCT, int nLenLinearCT, int *bOverflow )
{
    char szMult[16];
    int len1, len2;
    if ( mult > 0 && !*bOverflow && 0 < (len1 = strlen( szElement )) ) {
        if ( mult > 1 ) {
            len2 = sprintf( szMult, "%d", mult );
        } else {
            len2 = 0;
            szMult[0] = '\0';
        }
        if ( len1 + len2 < nLenLinearCT ) {
            memcpy( szLinearCT, szElement, len1 );
            memcpy( szLinearCT+len1, szMult, len2+1 ); /*  adding zero termination */
            return len1+len2;
        } else {
            (*bOverflow) ++;
        }
    }
    return 0;
}
/**********************************************************************************************/
/*  if num_C > 0 then nAtom does not contain C or H */
/*  otherwise all elements are in alphabetic order */
int MakeHillFormula( U_CHAR *nAtom, int num_atoms,
                  char *szLinearCT, int nLen_szLinearCT, int num_C, int num_H, int *bOverflow )
{
    char szElement[4];
    int  mult, compare2H;
    int  i, nLen, bOvfl;
    U_CHAR nPrevAtom;

    nLen       = 0;
    mult       = 0;
    bOvfl      = 0;
    nPrevAtom  = (U_CHAR)-2; /*  non-existent number */


    if ( num_C ) {
        nLen += AddElementAndCount( "C", num_C, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
        if ( num_H ) {
            nLen += AddElementAndCount( "H", num_H, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
            num_H = 0;
        }
    }

    for ( i = 0; i < num_atoms; i ++ ) {

        if ( nPrevAtom != nAtom[i] ) {
            if ( mult ) {
                nLen += AddElementAndCount( szElement, mult, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
            }
            mult = 1;
            if ( GetElementFormulaFromAtNum((int)nAtom[i], szElement ) ) {
                return -1; /*  wrong element */
            }
            nPrevAtom = nAtom[i];
            if ( !strcmp( "C", szElement ) ) {
                return -1;
            }
            compare2H = strcmp( "H", szElement );
            if ( !compare2H ) {
                return -1;
            }
            if ( compare2H < 0 && num_H ) {
                /*   H-atom should be located in front of szElement */
                nLen += AddElementAndCount( "H", num_H, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
                num_H = 0;
            }
        } else {
            mult ++;
        }
    }
    if ( mult ) {
        /*  the last element if any */
        nLen += AddElementAndCount( szElement, mult, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
    }
    if ( num_H ) {
        /*  if H has not been output... */
        nLen += AddElementAndCount( "H", num_H, szLinearCT+nLen, nLen_szLinearCT-nLen, &bOvfl );
    }
    *bOverflow |= (0 != bOvfl);
    return bOvfl? nLen_szLinearCT+1: nLen;
}
/**********************************************************************************************/
char *AllocateAndFillHillFormula( INChI *pINChI )
{
    int num_C, num_H, nLen, nNumNonHAtoms, ret, bOverflow;
    char *pHillFormula = NULL;
    bOverflow = 0;
    if ( !GetHillFormulaCounts( pINChI->nAtom, pINChI->nNum_H, pINChI->nNumberOfAtoms,
                          pINChI->nTautomer, pINChI->lenTautomer,
                          &num_C, &num_H, &nLen, &nNumNonHAtoms ) ) {
        if ( (pHillFormula = (char*) inchi_malloc( nLen+1 )) ) {
            ret = MakeHillFormula( pINChI->nAtom+num_C, nNumNonHAtoms-num_C,
                  pHillFormula, nLen+1, num_C, num_H, &bOverflow );
            if ( ret != nLen || bOverflow ) {
                inchi_free( pHillFormula );
                pHillFormula = NULL;
            }
        }
    }
    return pHillFormula;
}

/************************************************************************************/
/* return value: 0 => copied to stereo bonds; 1=> Allene copied to stereocenters    */
/* on input nNumberOfStereoBonds==NULL means second call, use Stereo->...Inv        */
/************************************************************************************/
int Copy2StereoBondOrAllene( INChI_Stereo *Stereo, int *nNumberOfStereoCenters, int *nNumberOfStereoBonds,
                            AT_STEREO_DBLE *LinearCTStereoDble,
                            AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic )
{
    int cumulene_len, j, next_j /* ordering number of the central allene atom */, next_neigh;
    AT_RANK  at_num;
    int      parity;
    if ( pCanonOrd && pCanonRank ) {
        j = pCanonOrd[(int)LinearCTStereoDble->at_num1-1];
        /* if allene then find the central atom, at[next_j] */
        if ( bIsotopic ) {
            cumulene_len = BOND_CHAIN_LEN(at[j].stereo_bond_parity2[0]);
            if ( cumulene_len % 2 && (1 >= MAX_NUM_STEREO_BONDS || !at[j].stereo_bond_neighbor2[1]) ) {
                next_j = at[j].neighbor[(int)at[j].stereo_bond_ord2[0]];
                for ( cumulene_len = (cumulene_len-1)/2; cumulene_len && 2==at[next_j].valence; cumulene_len -- ) {
                    next_neigh = (j == at[next_j].neighbor[0]);
                    j = next_j;
                    next_j = at[next_j].neighbor[next_neigh];
                }
                /* next_j is the central atom */
            } else {
                cumulene_len = -1; /* not an allene */
            }

        } else {
            cumulene_len = BOND_CHAIN_LEN(at[j].stereo_bond_parity[0]);
            if ( cumulene_len % 2 && (1 >= MAX_NUM_STEREO_BONDS || !at[j].stereo_bond_neighbor[1]) ) {
                next_j = at[j].neighbor[(int)at[j].stereo_bond_ord[0]];
                for ( cumulene_len = (cumulene_len-1)/2; cumulene_len && 2==at[next_j].valence; cumulene_len -- ) {
                    next_neigh = (j == at[next_j].neighbor[0]);
                    j = next_j;
                    next_j = at[next_j].neighbor[next_neigh];
                }
            } else {
                cumulene_len = -1; /* not an allene */
            }
        }
        if ( !cumulene_len ) {
            /* allene has been found; insert new stereocenter and parity */
            AT_NUMB *nNumber;
            S_CHAR  *t_parity;
            nNumber  = nNumberOfStereoBonds? Stereo->nNumber  : Stereo->nNumberInv;
            t_parity = nNumberOfStereoBonds? Stereo->t_parity : Stereo->t_parityInv;
            at_num = pCanonRank[next_j];
            parity = LinearCTStereoDble->parity;
            /* free room for the new stereocenter */
            for ( j = 0; j < *nNumberOfStereoCenters && Stereo->nNumber[j] < at_num; j ++ )
                ;
            if ( j < *nNumberOfStereoCenters ) {
                memmove( nNumber + j + 1, nNumber + j, (*nNumberOfStereoCenters-j)*sizeof(nNumber[0]) );
                memmove( t_parity + j + 1, t_parity + j, (*nNumberOfStereoCenters-j)*sizeof(t_parity[0]) );
            }
            /* fill the new stereo center info */
            
            nNumber[j] = at_num;
            t_parity[j] = parity;
            (*nNumberOfStereoCenters) ++;
            return 1;
        }
    }
    /* save the stereo bond info */
    if ( nNumberOfStereoBonds ) {
        j = *nNumberOfStereoBonds;
        Stereo->b_parity[j]   = LinearCTStereoDble->parity;
        Stereo->nBondAtom1[j] = LinearCTStereoDble->at_num1;
        Stereo->nBondAtom2[j] = LinearCTStereoDble->at_num2;
        (*nNumberOfStereoBonds) ++;
    }
    return 0;
}
/***************************************************************************/
int CopyLinearCTStereoToINChIStereo( INChI_Stereo *Stereo,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble
           , AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic
           , AT_STEREO_CARB *LinearCTStereoCarbInv
           , AT_STEREO_DBLE *LinearCTStereoDbleInv
           , AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv )
{
    int n, i, nErrorCode = 0, len;
    int bAllene;
    int diff;
    int lenInv, bAlleneInv;
    /*  stereo centers */
    n = Stereo->nNumberOfStereoCenters = nLenLinearCTStereoCarb;
    for ( i = 0; i < n; i ++ ) {
        Stereo->nNumber[i]  = LinearCTStereoCarb[i].at_num;
        Stereo->t_parity[i] = LinearCTStereoCarb[i].parity;
        Stereo->nNumberInv[i]  = LinearCTStereoCarbInv[i].at_num;
        Stereo->t_parityInv[i] = LinearCTStereoCarbInv[i].parity;
    }
    /*  stereo bonds */
    n = nLenLinearCTStereoDble;
    lenInv = Stereo->nNumberOfStereoCenters;
    for ( i = len = 0; i < n; i ++ ) {
        bAllene =
        Copy2StereoBondOrAllene( Stereo, &Stereo->nNumberOfStereoCenters,
                    &len, LinearCTStereoDble+i, pCanonOrd, pCanonRank, at, bIsotopic );
        bAlleneInv =
        Copy2StereoBondOrAllene( Stereo, &lenInv,
                    NULL, LinearCTStereoDbleInv+i, pCanonOrdInv, pCanonRankInv, at, bIsotopic );
        /* make sure double bond stereo is identical in original and inverted geometry */
        /* Note: all allenes are AFTER double bonds in LinearCTStereoDble... */
        if ( bAllene != bAlleneInv || (!bAllene &&
             CompareLinCtStereoDble ( LinearCTStereoDble+i,    1,
                                      LinearCTStereoDbleInv+i, 1 )) ) {
            nErrorCode = -4;          /* double bond stereo Inv is NOT identical to Abs */
            goto exit_function;
        }
    }
    Stereo->nNumberOfStereoBonds = len;

    if ( lenInv != Stereo->nNumberOfStereoCenters ) {
        nErrorCode = -5; /* different number of stereo centers in Abs and Inv */
        goto exit_function;
    }
    /* compare inverted stereocenters to absolute */
    n    = Stereo->nNumberOfStereoCenters;
    diff = 0;
    for ( i = 0, diff = 0; i < n; i ++ ) {
        if ( Stereo->nNumberInv[i] != Stereo->nNumber[i] ) {
            diff = (Stereo->nNumberInv[i] > Stereo->nNumber[i])? 2 : -2;
            break; /* Abs != Inv */
        }
        if ( Stereo->t_parityInv[i] != Stereo->t_parity[i] ) {
            diff = (Stereo->t_parityInv[i] > Stereo->t_parity[i])? 1 : -1;
            break; /* Abs != Inv */
        }
    }
    Stereo->nCompInv2Abs = (diff > 0)? 1 : (diff < 0)? -1 : 0;
    if ( diff == -1 || diff == 1 ) {
        /* the first found difference was in parities */
        for ( i = 0, diff = 0; i < n; i ++ ) {
            if ( Stereo->nNumberInv[i] != Stereo->nNumber[i] ) {
                diff = 2; /* difference in stereo center numbering */
                break;
            }
            /*  parities can be only 1, 2, 3, 4. Therefore only mutually inverted pairs
             *  (t_parityInv, t_parity) = (1,2) or (2,1) statisfy conditions
             *  (t_parityInv != t_parity) && (t_parityInv + t_parity == 3)
             */
            if ( Stereo->t_parityInv[i] == Stereo->t_parity[i] ||
                 Stereo->t_parityInv[i] +  Stereo->t_parity[i] != 3 ) {
                diff = 1; /* parities are same or different and cannot be obtained by simple inversion */
                break;
            }
        }
        Stereo->bTrivialInv = !diff;
    } else {
        Stereo->bTrivialInv = 0;
    }
exit_function:

    return nErrorCode;
}
/***************************************************************************/
int MarkAmbiguousStereo( sp_ATOM *at, inp_ATOM *norm_at, int bIsotopic, AT_NUMB *pCanonOrd,
           AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb,
           AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble )
{
    int n, i, j1, j2, num, mark_atom, mark_bond;
    
    if ( !pCanonOrd )
        return -1;
    num = 0;
    n = nLenLinearCTStereoCarb;
    mark_atom = bIsotopic? AMBIGUOUS_STEREO_ATOM_ISO : AMBIGUOUS_STEREO_ATOM;
    for ( i = 0; i < n; i ++ ) {
        /*  mark ambiguous stereo centers (for displaying and "Ambiguous stereo" message) */
        if ( ATOM_PARITY_NOT_UNKN(LinearCTStereoCarb[i].parity) &&
             at[j1=pCanonOrd[(int)LinearCTStereoCarb[i].at_num-1]].bAmbiguousStereo ) {
            at[j1].bAmbiguousStereo |= mark_atom;
            norm_at[j1].bAmbiguousStereo |= mark_atom;
            num ++;
        }
    }

    n = nLenLinearCTStereoDble;
    mark_bond = bIsotopic? AMBIGUOUS_STEREO_BOND_ISO : AMBIGUOUS_STEREO_BOND;
    for ( i = 0; i < n; i ++ ) {
        /*  mark ambiguous stereo bonds or allenes (for displaying and "Ambiguous stereo" message) */
        if ( ATOM_PARITY_WELL_DEF(LinearCTStereoDble[i].parity) ) {
            j1=pCanonOrd[(int)LinearCTStereoDble[i].at_num1-1];
            j2=pCanonOrd[(int)LinearCTStereoDble[i].at_num2-1];
            if ( at[j1].bAmbiguousStereo || at[j2].bAmbiguousStereo ) {
                /* if it is an allene then mark the central atom only
                   because the bonds should not be marked to avoid misleading
                   message "Ambiguous stereo: bond(s)": Allene makes a stereocenter
                */
                int j1_parity = bIsotopic? at[j1].stereo_bond_parity2[0] :
                                           at[j1].stereo_bond_parity[0];
                int cumulene_len = BOND_CHAIN_LEN(j1_parity); /* 0 => double bond, 1 => allene, 2 => cumulene,..*/
                if ( cumulene_len % 2 && (1 >= MAX_NUM_STEREO_BONDS ||
                     !(bIsotopic? at[j1].stereo_bond_neighbor2[1] :
                                  at[j1].stereo_bond_neighbor[1] )) ) {
                    /*  found an allene; locate its central atom */
                    int next_j, next_neigh;
                    int j = j1;
                    next_j =  at[j].neighbor[bIsotopic? at[j].stereo_bond_ord2[0] :
                                                        at[j].stereo_bond_ord[0] ];
                    for ( cumulene_len = (cumulene_len-1)/2;
                               cumulene_len && 2==at[next_j].valence;
                                      cumulene_len -- ) {
                        next_neigh = (j == at[next_j].neighbor[0]);
                        j = next_j;
                        next_j = at[next_j].neighbor[next_neigh];
                    }
                    /* next_j is the central atom */
                    if ( 2==at[next_j].valence ) {
                        at[next_j].bAmbiguousStereo |= mark_atom;
                        norm_at[next_j].bAmbiguousStereo |= mark_atom;
                        num ++;
                        continue; /* do not mark the cumulene "bond" endpoints */
                    }
                }
                /* not an allene, mark double bond or cumulene end atoms */
                if ( at[j1].bAmbiguousStereo ) {
                    at[j1].bAmbiguousStereo |= mark_bond; /*  ??? */
                    norm_at[j1].bAmbiguousStereo |= mark_bond;
                    num ++;
                }
                if ( at[j2].bAmbiguousStereo ) {
                    at[j2].bAmbiguousStereo |= mark_bond; /*  ??? */
                    norm_at[j2].bAmbiguousStereo |= mark_bond;
                    num ++;
                }
            }
        }
    }
    return num;

}
/**********************************************************************************************/
INCHI_MODE UnmarkAllUndefinedUnknownStereo( INChI_Stereo *Stereo, INCHI_MODE nUserMode )
{
    INCHI_MODE nRet = 0;
    int   i, n;
    if ( !Stereo || (Stereo && !Stereo->nNumberOfStereoCenters && !Stereo->nNumberOfStereoBonds)) {
        return nRet;
    }

    /* stereocenters */
    if ( !Stereo->nCompInv2Abs &&
         (n=Stereo->nNumberOfStereoCenters) > 0 && (nUserMode & REQ_MODE_SC_IGN_ALL_UU) ) {

        for ( i = 0; i < n && !ATOM_PARITY_WELL_DEF(Stereo->t_parity[i]); i ++ )
            ;
        if ( i == n ) {
            Stereo->nNumberOfStereoCenters = 0;
            for ( i = 0; i < n; i ++ ) {
                Stereo->t_parity[i] = 0;
                Stereo->nNumber[i]  = 0;
                Stereo->t_parityInv[i] = 0;
                Stereo->nNumberInv[i]  = 0;
            }
            nRet |= REQ_MODE_SC_IGN_ALL_UU;
        }
    }
    /* stereobonds */
    if ( (n=Stereo->nNumberOfStereoBonds) > 0 && (nUserMode & REQ_MODE_SB_IGN_ALL_UU) ) {
        for ( i = 0; i < n && !ATOM_PARITY_WELL_DEF(Stereo->b_parity[i]); i ++ )
            ;
        if ( i == n ) {
            Stereo->nNumberOfStereoBonds = 0;
            for ( i = 0; i < n; i ++ ) {
                Stereo->b_parity[i] = 0;
                Stereo->nBondAtom1[i]  = 0;
                Stereo->nBondAtom2[i]  = 0;
            }
            nRet |= REQ_MODE_SB_IGN_ALL_UU;
        }
    }

    return nRet;
}
#if ( defined(TARGET_API_LIB) || ADD_CMLPP==1 )
/**********************************************************************************************/
void WriteCoord( char *str, double x )
{
    if ( x < -9999999.9 ) {
        sprintf( str, "%10.2e", x );
    } else
    if ( x < -999999.99 ) {
        sprintf( str, "%10.2f", x );
    } else
    if ( x < -99999.999 ) {
        sprintf( str, "%10.3f", x );
    } else
    if ( x < 99999.9999 ) {
        sprintf( str, "%10.4f", x );
    } else
    if ( x < 999999.999 ) {
        sprintf( str, "%10.3f", x );
    } else
    if ( x < 9999999.99 ) {
        sprintf( str, "%10.2f", x );
    } else
    if ( x < 99999999.9 ) {
        sprintf( str, "%10.1f", x );
    } else {
        sprintf( str, "%10.3e", x );
    }
}
#endif
/* used CANON_STAT members

    pCS->LinearCT
    pCS->LinearCTIsotopic
    pCS->LinearCTIsotopicStereoCarb
    pCS->LinearCTIsotopicStereoCarbInv
    pCS->LinearCTIsotopicStereoDble
    pCS->LinearCTIsotopicStereoDbleInv
    pCS->LinearCTIsotopicTautomer
    pCS->LinearCTStereoCarb
    pCS->LinearCTStereoCarbInv
    pCS->LinearCTStereoDble
    pCS->LinearCTStereoDbleInv
    pCS->nCanonOrd
    pCS->nCanonOrdIsotopic
    pCS->nCanonOrdIsotopicStereo
    pCS->nCanonOrdIsotopicStereoInv
    pCS->nCanonOrdIsotopicStereoTaut
    pCS->nCanonOrdIsotopicTaut
    pCS->nCanonOrdStereo
    pCS->nCanonOrdStereoInv
    pCS->nCanonOrdStereoTaut
    pCS->nCanonOrdTaut
    pCS->nLenCanonOrd
    pCS->nLenCanonOrdIsotopic
    pCS->nLenCanonOrdIsotopicStereo
    pCS->nLenCanonOrdIsotopicStereoTaut
    pCS->nLenCanonOrdIsotopicTaut
    pCS->nLenCanonOrdStereo
    pCS->nLenCanonOrdStereoTaut
    pCS->nLenCanonOrdTaut
    pCS->nLenLinearCTAtOnly
    pCS->nLenLinearCTIsotopic
    pCS->nLenLinearCTIsotopicStereoCarb
    pCS->nLenLinearCTIsotopicStereoDble
    pCS->nLenLinearCTIsotopicTautomer
    pCS->nLenLinearCTStereoCarb
    pCS->nLenLinearCTStereoDble
    pCS->nNum_H
    pCS->nNum_H_fixed
    pCS->nSymmRank
    pCS->nSymmRankIsotopic
    pCS->nSymmRankIsotopicTaut
    pCS->nSymmRankTaut
    pCS->t_group_info
    pCS->t_group_info->num_t_groups

*/
/**********************************************************************************************/
int FillOutINChI( INChI *pINChI, INChI_Aux *pINChI_Aux,
                 int num_atoms, int num_at_tg, int num_removed_H,
                 sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, int bTautomeric,
                 INCHI_MODE nUserMode, char *pStrErrStruct )
{
    int i, j, m, n, g, len, ii, ret=0;

    AT_NUMB   *pSymmRank, *pOrigNosInCanonOrd, *pConstitEquNumb, *pCanonOrd=NULL, *pCanonOrdInv=NULL, *pCanonOrdTaut;
    T_GROUP_INFO     *t_group_info = pCS->t_group_info;
    T_GROUP *t_group;
    int nErrorCode = 0;
    AT_NUMB *pCanonRank, *pCanonRankInv; /* canonical ranks of the atoms or tautomeric groups */
    AT_NUMB *pCanonRankAtoms=NULL, *pSortOrd = NULL;
    AT_RANK nMinOrd;
    INChI_Stereo *Stereo;
    int          bUseNumberingInv = 0, bUseIsotopicNumberingInv = 0;
    INCHI_MODE    nStereoUnmarkMode;

    /*AT_NUMB  *pCanonOrdNonIso = NULL, *pCanonOrdIso = NULL;*/
    /*AT_NUMB  *nOrigAtNosInCanonOrdNonIso = NULL, *nOrigAtNosInCanonOrdIso = NULL;*/

    /*  Check for warnings */
    if ( pCS->nLenLinearCTStereoCarb < 0 || pCS->nLenLinearCTStereoDble  < 0 ||
         pCS->nLenCanonOrdStereo    < 0 || pCS->nLenCanonOrdStereoTaut < 0) {
        nErrorCode     |= WARN_FAILED_STEREO;
    }
    if ( pCS->nLenLinearCTIsotopic < 0  || pCS->nLenLinearCTIsotopicTautomer < 0 ||
         pCS->nLenCanonOrdIsotopic < 0 || pCS->nLenCanonOrdIsotopicTaut    < 0  ) {
        nErrorCode     |= WARN_FAILED_ISOTOPIC;
    }
    if ( pCS->nLenLinearCTIsotopicStereoCarb < 0 || pCS->nLenLinearCTIsotopicStereoDble  < 0 ||
         pCS->nLenCanonOrdIsotopicStereo    < 0 || pCS->nLenCanonOrdIsotopicStereoTaut < 0) {
        nErrorCode     |= WARN_FAILED_ISOTOPIC_STEREO;
    }
    pCanonRankAtoms = (AT_NUMB *)inchi_calloc( num_at_tg+1, sizeof(pCanonRankAtoms[0]) );
    pSortOrd        = (AT_NUMB *)inchi_calloc( num_at_tg+1, sizeof(pSortOrd[0]) ); /*  must have more than num_atoms */

    if ( !pCanonRankAtoms || !pSortOrd ) {
        nErrorCode = 0;
        ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = CT_OUT_OF_RAM;
        goto exit_function;
    }

    /*  total charge */
    for ( i = 0, n = 0; i < num_atoms+num_removed_H; i ++ ) {
        n += at[i].charge;
    }
    pINChI->nTotalCharge = n;

    /*  number of atoms */
    pINChI->nNumberOfAtoms     = num_atoms;
    pINChI_Aux->nNumberOfAtoms = num_atoms;

    /* removed protons and detachable isotopic H */
    if ( bTautomeric && t_group_info ) {
        pINChI_Aux->nNumRemovedProtons = t_group_info->tni.nNumRemovedProtons;
        for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
            pINChI_Aux->nNumRemovedIsotopicH[i] = t_group_info->num_iso_H[i] 
                                               + t_group_info->tni.nNumRemovedProtonsIsotopic[i];
        }
        if ( pINChI_Aux->bNormalizationFlags & FLAG_FORCE_SALT_TAUT ) {
            pINChI->nFlags |= INCHI_FLAG_HARD_ADD_REM_PROTON;
        }
        if ( pINChI_Aux->bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL) ) {
            AddMOLfileError(pStrErrStruct, "Proton(s) added/removed");
        }
        if ( pINChI_Aux->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL ) {
            AddMOLfileError(pStrErrStruct, "Charges neutralized");
        }
    }

    /* abs or rel stereo may establish one of two canonical numberings */
    if ( (pCS->nLenLinearCTStereoCarb > 0 || pCS->nLenLinearCTStereoDble > 0) &&
          pCS->nLenCanonOrdStereo > 0 &&
         ((pCS->LinearCTStereoCarb && pCS->LinearCTStereoCarbInv) ||
          (pCS->LinearCTStereoDble && pCS->LinearCTStereoDbleInv)) &&
          pCS->nCanonOrdStereo    && pCS->nCanonOrdStereoInv
       ) {

        pCanonRank    = pCanonRankAtoms;
        pCanonOrd     = pCS->nCanonOrdStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv  = pCS->nCanonOrdStereoInv;
        Stereo        = pINChI->Stereo;
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRankInv[pCanonOrdInv[i]] = 
            pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
                           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble
                           , pCanonOrd, pCanonRank, at, 0 /* non-isotopic */
                           , pCS->LinearCTStereoCarbInv
                           , pCS->LinearCTStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv ); 

        if ( Stereo->t_parityInv && Stereo->nNumberInv ) {
            if ( nUserMode & REQ_MODE_RELATIVE_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if ( nUserMode & REQ_MODE_RACEMIC_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if ( Stereo->nCompInv2Abs ) {
                if ( Stereo->nCompInv2Abs == -1 ) {
                    /* switch pointers in Stereo so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv  = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv  = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber     = nNumberInv;
                    Stereo->t_parity    = t_parityInv;
                    /* switch pointers to set rel. stereo to pINChI_Aux->nOrigAtNosInCanonOrd
                                       and inv. stereo to pINChI_Aux->nOrigAtNosInCanonOrdInv */
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
                    bUseNumberingInv    = 1; /* use inverted stereo numbering instead of normal */
                }
            }
        }

        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI_Aux->nOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nOrigAtNosInCanonOrd[i]    = at[pCanonOrd[i]].orig_at_number;
        }
        if ( bUseNumberingInv ) {
            /* switch ptrs back to avoid confusion */
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
            /* save inverted stereo ranks & order because it represents the smallest (relative) */
            memcpy( pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]) );
            /* change pCS->nCanonOrdStereo[] to inverted: */
            memcpy( pCanonOrd,  pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]) );
        }
        pCanonRankInv = NULL;
        pCanonOrdInv  = NULL;
        pOrigNosInCanonOrd = NULL;

    } else { /*------------------------------ no stereo */

        pCanonOrd            = pCS->nLenCanonOrdStereo > 0? pCS->nCanonOrdStereo :
                               pCS->nLenCanonOrd       > 0? pCS->nCanonOrd : NULL;
        pCanonRank           = pCanonRankAtoms;
        pOrigNosInCanonOrd   = pINChI_Aux->nOrigAtNosInCanonOrd;
        if ( pCanonOrd && pCanonRank ) {
            for ( i = 0; i < num_atoms; i ++ ) {
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
                pOrigNosInCanonOrd[i]          = at[pCanonOrd[i]].orig_at_number;
            }
            for ( ; i < num_at_tg; i ++ ) {
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
            }
        }
    }
    /*pCanonOrdNonIso = pCanonOrd;*/  /* save for aux info */


    if ( pINChI_Aux->OrigInfo ) {
        /* charges, radicals, valences */
        for ( i = 0; i < num_atoms; i ++ ) {
            ii = pCanonOrd[i];
            if ( norm_at[ii].valence || norm_at[ii].num_H ) {
                pINChI_Aux->OrigInfo[i].cCharge  = norm_at[ii].charge;
                pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical==RADICAL_SINGLET)? 0 :
                                                  (norm_at[ii].radical==RADICAL_DOUBLET)? 1 :
                                                  (norm_at[ii].radical==RADICAL_TRIPLET)? 2 :
                                                  norm_at[ii].radical? 3 : 0 ;
                pINChI_Aux->OrigInfo[i].cUnusualValence = 
                    get_unusual_el_valence( norm_at[ii].el_number, norm_at[ii].charge, norm_at[ii].radical,
                                            norm_at[ii].chem_bonds_valence, norm_at[ii].num_H, norm_at[ii].valence );
            } else {
                /* charge of a single atom component is in the INChI; valence = 0 is standard */
                pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical==RADICAL_SINGLET)? 0 :
                                                  (norm_at[ii].radical==RADICAL_DOUBLET)? 1 :
                                                  (norm_at[ii].radical==RADICAL_TRIPLET)? 2 :
                                                  norm_at[ii].radical? 3 : 0 ;
            }

        }
    }

    /* non-isotopic canonical numbers and equivalence of atoms (Aux) */
    pConstitEquNumb      = pINChI_Aux->nConstitEquNumbers;  /*  contitutional equivalence */
    pSymmRank            = pCS->nSymmRank;
    if ( pCanonOrd && pCanonRank && pSymmRank && pConstitEquNumb ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pConstitEquNumb[i]       = pSymmRank[pCanonOrd[i]]; /*  constit. equ. ranks in order of canonical numbers */
            pSortOrd[i]              = i;
        }
        for ( ; i < num_at_tg; i ++ ) {
            pSortOrd[i]              = MAX_ATOMS; /* for debugging only */
        }
        pn_RankForSort  = pConstitEquNumb;
        qsort( pSortOrd, num_atoms, sizeof(pSortOrd[0]), CompRanksOrd );
        for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j ++ ) {
            if ( j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                nMinOrd ++;
                if ( j - i > 1 ) {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while ( i < j ) {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd; /*  = min. canon. rank in the group of equ. atoms */
                    }
                    /*  at this point j == i */
                } else {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /*  means the atom is not equivalent to any other */
                }
                nMinOrd = pSortOrd[j]; /*  at the end j = num_atoms */
            }
        }
    } else {
        nErrorCode  |= ERR_NO_CANON_RESULTS;
        ret = -1;  /*  program error; no breakpoint here */
        goto exit_function;
    }
    /*  atomic numbers from the Periodic Table */
    for ( i = 0; i < num_atoms; i ++ ) {
        pINChI->nAtom[i] = (int)at[pCanonOrd[i]].el_number;
    }
    /*  connection table: atoms only (before 7-29-2003 pCS->LinearCT2 contained non-isotopic CT) */
    if ( pCS->nLenLinearCTAtOnly <= 0 || !pCS->LinearCT || !pINChI->nConnTable ) {
        nErrorCode  |= ERR_NO_CANON_RESULTS;
        ret = -2;
        goto exit_function;
    }
    memcpy( pINChI->nConnTable, pCS->LinearCT, sizeof(pINChI->nConnTable[0])*pCS->nLenLinearCTAtOnly);
    pINChI->lenConnTable = pCS->nLenLinearCTAtOnly;
    
    /*  tautomeric group(s) canonical representation */
    len = 0;
    if ( bTautomeric && 0 < (n = SortTautomerGroupsAndEndpoints( t_group_info, num_atoms, num_at_tg, pCanonRank )) ) {
        /* SortTautomerGroupsAndEndpoints() produces canonically ordered t-groups */
        pINChI->nFlags |= (t_group_info->bTautFlagsDone & TG_FLAG_ALL_SALT_DONE)? INCHI_FLAG_ACID_TAUT : 0;
        /*  number of tautomeric groups */
        pINChI->nTautomer[len ++] = (AT_NUMB)n;
        /* store each tautomeric group, one by one */
        for ( i = 0; i < n; i ++ ) {
            g = (int)t_group_info->tGroupNumber[i]; /* original group numbers in sorted order */
            t_group = t_group_info->t_group + g;    /* pointer to the tautomeric group */
            /*  NumAt+INCHI_T_NUM_MOVABLE (group length excluding this number) */
            pINChI->nTautomer[len ++]     = t_group->nNumEndpoints+INCHI_T_NUM_MOVABLE;
            /*  Num(H), Num(-) */
            for ( j = 0; j < INCHI_T_NUM_MOVABLE && j < T_NUM_NO_ISOTOPIC; j ++ )
                pINChI->nTautomer[len ++]     = t_group->num[j];
            for ( j = T_NUM_NO_ISOTOPIC; j < INCHI_T_NUM_MOVABLE; j ++ )
                pINChI->nTautomer[len ++]     = 0; /* should not happen */
            /* tautomeric group endpoint canonical numbers, pre-sorted in ascending order */
            for ( j  = (int)t_group->nFirstEndpointAtNoPos,
                  m  = j + (int)t_group->nNumEndpoints; j < m; j ++ ) {
                pINChI->nTautomer[len ++] = pCanonRank[(int)t_group_info->nEndpointAtomNumber[j]]; /*  At[j] */
            }
        }
        pINChI->lenTautomer = len;
        pINChI_Aux->nNumberOfTGroups = n;
    } else {
        pINChI->lenTautomer = 0;
        pINChI_Aux->nNumberOfTGroups = 0;
        if ( t_group_info && ((t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) ||
                              (t_group_info->nNumIsotopicEndpoints>1 &&
                              (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))))
           ) {
            /* only protons (re)moved or added */
            pINChI->lenTautomer  = 1;
            pINChI->nTautomer[0] = 0;
        }
    }

    /*  number of H (excluding tautomeric) */
    if ( pCS->nNum_H ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI->nNum_H[i] = pCS->nNum_H[i];
        } 
    }
    /*  number of fixed H (tautomeric H in non-tautomeric representation) */
    if ( pCS->nNum_H_fixed && !pINChI->lenTautomer ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI->nNum_H_fixed[i]  = pCS->nNum_H_fixed[i];
            pINChI->nNum_H[i]       += pCS->nNum_H_fixed[i];
        } 
    }

    /***********************************************************
     *  tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if ( pINChI->lenTautomer && (n=pINChI_Aux->nNumberOfTGroups) ) {
        pCanonOrdTaut   = pCS->nLenCanonOrdStereoTaut > 0? pCS->nCanonOrdStereoTaut :
                          pCS->nLenCanonOrdTaut       > 0? pCS->nCanonOrdTaut : NULL;
        pConstitEquNumb = pINChI_Aux->nConstitEquTGroupNumbers;
        pSymmRank       = pCS->nSymmRankTaut;
        if ( pCanonOrdTaut && pSymmRank && pConstitEquNumb ) {
            for ( i = 0; i < n; i ++ ) {
                pConstitEquNumb[i]       = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i]              = i;
            }
            pn_RankForSort  = pConstitEquNumb;
            qsort( pSortOrd, n, sizeof(pSortOrd[0]), CompRanksOrd );
            for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j ++ ) {
                if ( j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                    nMinOrd ++; /* make is start from 1, not from zero */
                    if ( j - i > 1 ) {
                        /*  found a sequence of more than one equivalent t-groups: i..j-1 */
                        while ( i < j ) {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    } else {
                        pConstitEquNumb[pSortOrd[i++]] = 0;
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j == n */
                }
            }
        }
    }

    /*  Allocate and fill Hill formula */
    if ( !(pINChI->szHillFormula = AllocateAndFillHillFormula( pINChI ) ) ) {
        nErrorCode = 0;
        ret = CT_WRONG_FORMULA; /* CT_OUT_OF_RAM;*/  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = ret;
        goto exit_function;
    }

    if ( (nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->Stereo, nUserMode )) ) {
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_UU : 0;    
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU)? INCHI_FLAG_SB_IGN_ALL_UU : 0;
        if ( (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU) ||
             (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU) ) {
             AddMOLfileError(pStrErrStruct, "Omitted undefined stereo"); 
        }
    }

    /*************************/
    /* mark ambiguous stereo */
    /*************************/
    MarkAmbiguousStereo( at, norm_at, 0 /* non-isotopic */, pCanonOrd,
           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble );


    /************************************************************************
     *
     *  isotopic part
     */
    /* abs or rel stereo may establish one of two canonical numberings */
    if ( (pCS->nLenLinearCTIsotopicStereoCarb > 0 || pCS->nLenLinearCTIsotopicStereoDble > 0) &&
          pCS->nLenCanonOrdIsotopicStereo > 0 &&
         ((pCS->LinearCTIsotopicStereoCarb && pCS->LinearCTIsotopicStereoCarbInv) ||
          (pCS->LinearCTIsotopicStereoDble && pCS->LinearCTIsotopicStereoDbleInv)) &&
          pCS->nCanonOrdIsotopicStereo    && pCS->nCanonOrdIsotopicStereoInv
          ) {
        /* found isotopic stereo */
        pCanonRank    = pCanonRankAtoms;
        pCanonOrd     = pCS->nCanonOrdIsotopicStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv  = pCS->nCanonOrdIsotopicStereoInv;
        Stereo        = pINChI->StereoIsotopic;
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRankInv[pCanonOrdInv[i]] =
            pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
                           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble
                           , pCanonOrd, pCanonRank, at, 1 /* isotopic */
                           , pCS->LinearCTIsotopicStereoCarbInv
                           , pCS->LinearCTIsotopicStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv ); 

        if ( Stereo->t_parityInv && Stereo->nNumberInv ) {
            if ( nUserMode & REQ_MODE_RELATIVE_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if ( nUserMode & REQ_MODE_RACEMIC_STEREO ) {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if ( Stereo->nCompInv2Abs ) {
                if ( Stereo->nCompInv2Abs == -1 ) {
                    /* switch pointers so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv  = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv  = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber     = nNumberInv;
                    Stereo->t_parity    = t_parityInv;
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
                    bUseIsotopicNumberingInv    = 1;
                }
            }
        }

        for ( i = 0; i < num_atoms; i ++ ) {
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[i]    = at[pCanonOrd[i]].orig_at_number;
        }
        if ( bUseIsotopicNumberingInv ) {
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd,  &pCanonOrdInv  );
            memcpy( pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]) );
            memcpy( pCanonOrd,  pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]) );
        }
        pCanonRankInv = NULL;
        pCanonOrdInv  = NULL;
        pOrigNosInCanonOrd = NULL;

    } else {
        /* no isotopic stereo */
        pCanonOrd = pCS->nLenCanonOrdIsotopicStereo > 0? pCS->nCanonOrdIsotopicStereo :
                    pCS->nLenCanonOrdIsotopic       > 0? pCS->nCanonOrdIsotopic : NULL;
        pCanonRank           = pCanonRankAtoms;
        pOrigNosInCanonOrd   = pINChI_Aux->nIsotopicOrigAtNosInCanonOrd;
        if ( pCanonOrd && pCanonRank ) {
            for ( i = 0; i < num_atoms; i ++ ) { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
                pOrigNosInCanonOrd[i]          = at[pCanonOrd[i]].orig_at_number;
            }
            for ( ; i < num_at_tg; i ++ ) { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]]       = (AT_NUMB)(i+1);
            }
        }
    }
    /*pCanonOrdIso = pCanonOrd;*/

    pConstitEquNumb      = pINChI_Aux->nConstitEquIsotopicNumbers;
    pSymmRank            = pCS->nSymmRankIsotopic;
    if ( pCanonOrd && pCanonRank && pConstitEquNumb && pSymmRank ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            pConstitEquNumb[i]       = pSymmRank[pCanonOrd[i]];
            pSortOrd[i]              = i;
        }
        for ( ; i < num_at_tg; i ++ ) {
            pSortOrd[i]              = i;
        }
        pn_RankForSort  = pConstitEquNumb;
        qsort( pSortOrd, num_atoms, sizeof(pSortOrd[0]), CompRanksOrd );
        for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j ++ ) {
            if ( j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                nMinOrd ++;
                if ( j - i > 1 ) {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while ( i < j ) {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                    }
                } else {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /* nMinOrd; */
                }
                nMinOrd = pSortOrd[j];
            }
        }
    } else {
        goto exit_function; /*  no isotopic info available */
    }
    /*  isotopic atoms */
    n = pINChI->nNumberOfIsotopicAtoms = pCS->nLenLinearCTIsotopic;
    for ( i = 0; i < n; i ++ ) {
        pINChI->IsotopicAtom[i].nAtomNumber    = pCS->LinearCTIsotopic[i].at_num;
        pINChI->IsotopicAtom[i].nIsoDifference = pCS->LinearCTIsotopic[i].iso_atw_diff;
        pINChI->IsotopicAtom[i].nNum_H         = pCS->LinearCTIsotopic[i].num_1H;
        pINChI->IsotopicAtom[i].nNum_D         = pCS->LinearCTIsotopic[i].num_D;
        pINChI->IsotopicAtom[i].nNum_T         = pCS->LinearCTIsotopic[i].num_T;
    }
    /*  isotopic tautomeric groups */
    n = pINChI->nNumberOfIsotopicTGroups = pCS->nLenLinearCTIsotopicTautomer;
    for ( i = 0; i < n; i ++ ) {
        pINChI->IsotopicTGroup[i].nTGroupNumber = pCS->LinearCTIsotopicTautomer[i].tgroup_num;
        pINChI->IsotopicTGroup[i].nNum_H        = pCS->LinearCTIsotopicTautomer[i].num[2];
        pINChI->IsotopicTGroup[i].nNum_D        = pCS->LinearCTIsotopicTautomer[i].num[1];
        pINChI->IsotopicTGroup[i].nNum_T        = pCS->LinearCTIsotopicTautomer[i].num[0];
    }
    /* atoms that may exchange isotopic H-atoms */
    if ( pCS->nExchgIsoH && pINChI->nPossibleLocationsOfIsotopicH ) {
        for ( i = 0, j = 1; i < num_atoms; i ++ ) {
            if ( pCS->nExchgIsoH[i] ) {
                pINChI->nPossibleLocationsOfIsotopicH[j++] = (AT_NUMB)(i+1); /* canonical number */
            }
        }
        pINChI->nPossibleLocationsOfIsotopicH[0] = (AT_NUMB)j; /* length including the 0th element */
    }

    if ( (nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->StereoIsotopic, nUserMode )) ) {
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;    
        pINChI->nFlags |= (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU)? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;    
        if ( (nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU) ||
             (nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU) ) {
             AddMOLfileError(pStrErrStruct, "Omitted undefined stereo"); 
        }
    }
    /* mark ambiguous stereo */
    MarkAmbiguousStereo( at, norm_at, 1 /* isotopic */, pCanonOrd,
           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble );

    /***********************************************************
     *  isotopic tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if ( pINChI->lenTautomer && pINChI_Aux->nConstitEquIsotopicTGroupNumbers && pCS->nSymmRankIsotopicTaut &&
         (pCS->nLenLinearCTIsotopic || pCS->nLenLinearCTIsotopicTautomer) &&
          t_group_info && t_group_info->num_t_groups > 0 ) {
        n = t_group_info->num_t_groups;
        pCanonOrdTaut   = pCS->nLenCanonOrdIsotopicStereoTaut > 0? 
                              (n=pCS->nLenCanonOrdIsotopicStereoTaut, pCS->nCanonOrdIsotopicStereoTaut) :
                          pCS->nLenCanonOrdIsotopicTaut       > 0?
                              (n=pCS->nLenCanonOrdIsotopicTaut,pCS->nCanonOrdIsotopicTaut) : (n=0,(AT_RANK*)NULL);
        pConstitEquNumb = pINChI_Aux->nConstitEquIsotopicTGroupNumbers;
        pSymmRank       = pCS->nSymmRankIsotopicTaut;
        if ( pCanonOrdTaut && pSymmRank && pConstitEquNumb && n > 0 ) {
            for ( i = 0; i < n; i ++ ) {
                pConstitEquNumb[i]       = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i]              = i;
            }
            pn_RankForSort  = pConstitEquNumb;
            qsort( pSortOrd, n, sizeof(pSortOrd[0]), CompRanksOrd );
            for ( i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j ++ ) {
                if ( j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]] ) {
                    nMinOrd ++;
                    if ( j - i > 1 ) {
                        /*  found a sequence of equivalent t-groups: i..j-1 */
                        while ( i < j ) {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    } else {
                        pConstitEquNumb[pSortOrd[i++]] = 0; /*  nMinOrd; */
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j = n */
                }
            }
        }
    }


exit_function:
    if ( pCanonRankAtoms )
        inchi_free( pCanonRankAtoms );
    if ( pSortOrd )
        inchi_free( pSortOrd );

    pINChI->nErrorCode     |= nErrorCode;
    pINChI_Aux->nErrorCode |= nErrorCode;

    return ret;
}
