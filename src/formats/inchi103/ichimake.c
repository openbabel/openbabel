/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#include "mode.h"

#if( TEST_RENUMB_ATOMS == 1 )
#include "ichitime.h"
#endif
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
#include "ichister.h"
#include "ichi_io.h"

int inp2spATOM( inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at );
int GetElementAndCount( const char **f, char *szEl, int *count );
int CompareHillFormulas( const char *f1, const char *f2 );
int CompareInchiStereo( INChI_Stereo *Stereo1, INCHI_MODE nFlags1, INChI_Stereo *Stereo2, INCHI_MODE nFlags2 );
int CompareReversedStereoINChI( INChI_Stereo *s1/* InChI from reversed struct */, INChI_Stereo *s2 /* input InChI */);
int GetAtomOrdNbrInCanonOrd( inp_ATOM *norm_at, AT_NUMB *nAtomOrdNbr,
                            AT_NUMB *nOrigAtNosInCanonOrd, int num_at );
int FillOutCanonInfAtom(inp_ATOM *norm_at, INF_ATOM_DATA *inf_norm_at_data, int init_num_at, int bIsotopic,
                        INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode);
int FillOutOneCanonInfAtom(inp_ATOM *inp_norm_at, INF_ATOM_DATA *inf_norm_at_data,
                           AT_NUMB *pStereoFlags, int init_num_at, int offset, int offset_H, int bIsotopic,
                           INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode);
int FillOutInputInfAtom(inp_ATOM *inp_at, INF_ATOM_DATA *inf_at_data, int init_num_at, int num_removed_H,
                        int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, int bAbcNumbers);

int CheckCanonNumberingCorrectness(
                 int num_atoms, int num_at_tg,
                 sp_ATOM *at, CANON_STAT *pCS, int bTautomeric,
                 char *pStrErrStruct );

static int CompareDfsDescendants4CT( const void *a1, const void *a2 );
int GetSp3RelRacAbs( const INChI *pINChI, INChI_Stereo *Stereo );


#if( TEST_RENUMB_ATOMS == 1 || READ_INCHI_STRING == 1 ) /*  { */
int CompareStereoINChI( INChI_Stereo *s1, INChI_Stereo *s2 );
#endif

#if( READ_INCHI_STRING == 1 ) /*  { */
/*************************************************************************************/

int CompareReversedStereoINChI2( INChI_Stereo *s1, INChI_Stereo *s2, ICR *picr);

#endif
/**********************************************************************************************/
int inp2spATOM( inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at )
{
    int i, j, val;
    memset( at, 0, sizeof(at[0])*num_inp_at );
    for ( i = 0; i < num_inp_at; i ++ ) {
        strncpy( at[i].elname, inp_at[i].elname, sizeof(at[0].elname) );
        at[i].el_number = (U_CHAR)get_periodic_table_number( at[i].elname );
        val = at[i].valence = inp_at[i].valence;
        for ( j = 0; j < val; j ++ ) {
            at[i].neighbor[j]  = inp_at[i].neighbor[j];
            at[i].bond_type[j] = inp_at[i].bond_type[j];
        }
        at[i].chem_bonds_valence = inp_at[i].chem_bonds_valence;
        at[i].orig_at_number     = inp_at[i].orig_at_number;
        at[i].orig_compt_at_numb= inp_at[i].orig_compt_at_numb;
        at[i].endpoint           = inp_at[i].endpoint;
        at[i].iso_atw_diff       = inp_at[i].iso_atw_diff;
        at[i].num_H              = inp_at[i].num_H;
        at[i].cFlags             = inp_at[i].cFlags;
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
            at[i].num_iso_H[j]       = inp_at[i].num_iso_H[j];
        }
        at[i].charge             = inp_at[i].charge;
        at[i].radical            = inp_at[i].radical;

#if( FIND_RING_SYSTEMS == 1 )
        at[i].nBlockSystem       = inp_at[i].nBlockSystem;
        at[i].bCutVertex         = inp_at[i].bCutVertex;
        at[i].nRingSystem        = inp_at[i].nRingSystem;
        at[i].nNumAtInRingSystem = inp_at[i].nNumAtInRingSystem;
#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
        at[i].nDistanceFromTerminal = inp_at[i].nDistanceFromTerminal;
#endif
#endif

/*
        at[i].x                  = inp_at[i].x;
        at[i].y                  = inp_at[i].y;
        at[i].z                  = inp_at[i].z;
*/
    }
    return 0;
}
/**********************************************************************************************/
int GetElementAndCount( const char **f, char *szEl, int *count )
{
    const char *p = *f;
    char       *q;
    int   i = 0;
    if ( *p ) {
        if ( isupper( UCINT *p ) ) {
            szEl[i++] = *p++;
            if ( *p && islower( UCINT *p ) ) {
                szEl[i++] = *p++;
            }
            szEl[i] = '\0';
            if ( 1 == i && szEl[0] == 'C' ) {
                szEl[0] = 'A'; /*  less than any element: */
                               /*  carbon-containing compounds should be first */
            }
            if ( *p && isdigit( UCINT *p ) ) {
                *count = strtol( p, &q, 10 );
                p = q;
            } else {
                *count = 1;
            }
            *f = p; /*  next element; */
            return 1;
        }
        return -1; /*  not a chemical formula */
    }
    strcpy( szEl, "Zz" ); /*  zero termination 'element' is larger than any other element */
    *count = 9999;        /*  zero termination 'element count' is larger than any other count */
    return 0;
}
/**********************************************************************************************

   E1 < E2 if strcmp( E1, E2) < 0  OR E2 is empty and E1 is not
   n1 < n2 if value n1 > n2

   Sorting order:

   C10H22N
   C10H22
   C2
   Ag2Cl2
   Ag2Cl
   Ag2F2
   Ag2
   AgCl
   AgF
   F6S
   F2S

**********************************************************************************************/
int CompareHillFormulas( const char *f1, const char *f2 )
{
    char szEl1[4], szEl2[4];
    int  count1, count2, ret1, ret2, ret;

    do {
        ret1 = GetElementAndCount( &f1, szEl1, &count1 );
        ret2 = GetElementAndCount( &f2, szEl2, &count2 );
        if ( 0 <= ret1 && 0 <= ret2 ) {
            if ( ret = strcmp( szEl1, szEl2 ) ) {
                return ret; /*  lexicographic order, string termination > any character */
            }
            if ( ret = count2 - count1 ) {
                return ret; /*  inverse atom count order */
            }
        } else {
            return 0; /*  program error <BRKPT> */
        }

    } while ( 0 < ret1 && 0 < ret2 );

    return 0;
}
/**********************************************************************************************/
int CompareHillFormulasNoH( const char *f1, const char *f2, int *num_H1, int *num_H2 )
{
    char szEl1[4], szEl2[4];
    int  count1, count2, ret1, ret2, ret;

    do {
        ret1 = GetElementAndCount( &f1, szEl1, &count1 );
        if ( 0 < ret1 && szEl1[0] == 'H' && !szEl1[1] ) {
            *num_H1 += count1;
            ret1 = GetElementAndCount( &f1, szEl1, &count1 );
        }
        ret2 = GetElementAndCount( &f2, szEl2, &count2 );
        if ( 0 < ret2 && szEl2[0] == 'H' && !szEl2[1] ) {
            *num_H2 += count2;
            ret2 = GetElementAndCount( &f2, szEl2, &count2 );
        }
        if ( 0 <= ret1 && 0 <= ret2 ) {
            if ( ret = strcmp( szEl1, szEl2 ) ) {
                return ret; /*  lexicographic order, string termination > any character */
            }
            if ( ret = count2 - count1 ) {
                return ret; /*  inverse atom count order */
            }
        } else {
            return 0; /*  program error <BRKPT> */
        }

    } while ( 0 < ret1 && 0 < ret2 );

    return 0;
}

/**************************************************************/
int CompareTautNonIsoPartOfINChI( const INChI *i1, const INChI *i2 )
{
    int len1, len2, ret, i;

    len1 = i1->lenTautomer > 0 && i1->nTautomer[0]? i1->lenTautomer:0;
    len2 = i2->lenTautomer > 0 && i2->nTautomer[0]? i2->lenTautomer:0;
    if ( ret  = len2 - len1 ) {
        return ret;
    }
    for ( i = 0; i < len1; i ++ ) {
        if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
            return ret;
    }
    return 0;
}

/**********************************************************************************************/
/*  sorting in descending order: return -1 if *p1 > *p2, return +1 if *p1 < *p2               */
/**********************************************************************************************/
int CompINChITautVsNonTaut(const INCHI_SORT *p1, const INCHI_SORT *p2, int bCompareIsotopic)
{
    int ret, num, i, num_H1, num_H2;
    
    const INChI *i1  = NULL; /* Mobile-H layers in Mobile-H sorting order */
    const INChI *i2  = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    
    int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    
    /* INChI_Stereo *Stereo1, *Stereo2; */

    n1 = ( p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms )? TAUT_YES : TAUT_NON;

    i1  = p1->pINChI[n1];
    i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
           p2->pINChI[TAUT_NON]->nNumberOfAtoms)? p2->pINChI[TAUT_NON] : (const INChI *)NULL;


    /* non-deleted-non-empty < deleted < empty */
    if ( i1 && !i2 )
        return 0;   /* non-empty is the smallest (first) */
    if ( !i1 && i2 )
        return 0;
    if ( !i1 && !i2 )
        return 0;
    if ( i1->bDeleted )
        return 1;    /* deleted is the largest (last) among non-empty */
    if ( i2->bDeleted )
        return -1;

    if ( i1->nNumberOfAtoms > 0 && !i2->nNumberOfAtoms )
        return 0;

    i2 = i2;

    num_H1 = num_H2 = 0;
    
    /* do not compare terminal H */
    if ( ret = CompareHillFormulasNoH( i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2 ) ) {
        return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    }

    /*********************************************************
            compare non-isotopic non-tautomeric part
     *********************************************************/

    /* compare number of atoms (excluding terminal H) */
    if ( ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms )
        return ret; /*  more atoms first */
    
    /*  compare elements  (excluding terminal H) */
    num = i1->nNumberOfAtoms;
    for ( i = 0; i < num; i ++ ) { /* should always be equal if Hill formulas are same */
        if ( ret = (int)i2->nAtom[i] - (int)i1->nAtom[i] )
            return ret; /* greater periodic number first */
    }
    /**********************************************************
        compare connection tables
    ***********************************************************/
    if ( ret = i2->lenConnTable - i1->lenConnTable )
        return ret; /* longer connection table first */
    num = i2->lenConnTable; 
    for ( i = 0; i < num; i ++ ) {
        if ( ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i] )
            return ret; /* greater connection table first */
    }
    /*********************************************************
      compare compare total number of H (inverse: H3 < H2 )
    **********************************************************/
    if ( ret = num_H2 - num_H1 )
        return ret;
    /*********************************************************
      compare non-tautomeric num_H: N < NH3 < NH2 < NH
    **********************************************************/
    num = i1->nNumberOfAtoms;
    for ( i = 0; i < num; i ++ ) {
        if ( i2->nNum_H[i] != i1->nNum_H[i] ) {
            return !i2->nNum_H[i]?  1 :  /* no H first */
                   !i1->nNum_H[i]? -1 :
                   (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
        }
    }
    /*********************************************************
         compare non-isotopic tautomeric part
     *********************************************************/
    if ( ret = CompareTautNonIsoPartOfINChI( i1, i2) ) {
        return ret;
    }
    /*
    if ( ret = i2->lenTautomer - i1->lenTautomer )
        return ret;
    num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    for ( i = 0; i < num; i ++ ) {
        if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
            return ret;
    }
    */
    /*********************************************************
     *                                                       *
     *  at this point both components are either tautomeric  *
     *  or non-tautomeric                                    *
     *                                                       *
     *********************************************************/

    /*********************************************************
       non-tautomeric "fixed H" specific
     *********************************************************/
    if ( /*TAUT_NON == bTaut &&*/ (i2 && i2->nNum_H_fixed ) ) {
        /* first, compare non-tautomeric chem. formulas -- they may be different */
        /* secondly, compare fixed-H distribution */
        if ( i2->nNum_H_fixed ) {
            num = i2->nNumberOfAtoms;
            for ( i = 0; i < num; i ++ ) {
                if ( i2->nNum_H_fixed[i] != 0 ) {
                    return 1;
                }
            }
        }
    }
    /*********************************************************
        compare non-isotopic stereo
     *********************************************************/
    ret = CompareInchiStereo( i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags );
    if ( ret ) {
        return ret;
    }
    /*******************************************************
        do not switch back to tautomeric i1, i2
     *******************************************************/
    /* -- how to switch back --
    if ( i1t ) {
        i1  = i1t;
        i1t = NULL;
    }
    if ( i2t ) {
        i2  = i2t;
        i2t = NULL;
    }
    */
    /******************************************************
         compare isotopic non-tautomeric part
     ******************************************************/
    if ( bCompareIsotopic ) {
        if ( ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms )
            return ret;
        num = i1->nNumberOfIsotopicAtoms;
        /*  compare isotopic atoms */
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference )
                return ret;
        }
        /* compare isotopic H */
        /* if tautomeric comparison mode then here are compared only non-tautomeric H */
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H )
                return ret;
        }
        /*****************************************************
             compare isotopic tautomeric part
         *****************************************************/
        if ( ret = i2->nNumberOfIsotopicTGroups || i1->nNumberOfIsotopicTGroups )
            return ret;
        /*
        num = i1->nNumberOfIsotopicTGroups;
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H )
                return ret;
        }
        */
        /****************************************************
            compare isotopic stereo
         ****************************************************/
        ret = CompareInchiStereo( i1->StereoIsotopic, i1->nFlags, i2->StereoIsotopic, i2->nFlags );
        if ( ret ) {
            return ret;
        }

    }


    /**********************************************************
        compare charges: non-charged first, then in order of
        ascending charges (negative first) 
    ***********************************************************/
    if ( i2->nTotalCharge && i1->nTotalCharge ) {
        /*  both are charged; smaller charges first */
        ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
        return ret;
    }
    if ( ret = (i1->nTotalCharge? 1:0) - (i2->nTotalCharge? 1:0) ) {
        /*  only one is charged; uncharged first */
        return ret;
    }
    /* stable sort */
    /*ret = p1->ord_number - p2->ord_number;*/

    return ret;
}

/*************************** stereo ***********************************************************/
typedef enum tagSp3StereoTypeTmp {
    SP3_NONE = 0,  /* no sp3 stereo: no /t, /m, /s segments */
    /* /t is present: */
    SP3_ONLY = 1,  /* no /s or /m segment: inversion leaves the structure unchanged */
    SP3_ABS  = 2,  /* abs stereo: both /m and /s are present */
    SP3_REL  = 4,  /* rel stereo: /s is present, /m is not */
    SP3_RAC  = 8,  /* racemic stereo: /s is presen, /m is nott */
    SP3_TYPE = (SP3_ABS|SP3_REL|SP3_RAC),          /* bitmap for checking the presence of /m */
    SP3_ANY  = (SP3_ABS|SP3_REL|SP3_RAC|SP3_ONLY)  /* bitmap for checking the presence of /t */
} SP3_TYPE_TMP;

/**********************************************************************************************/
int GetSp3RelRacAbs( const INChI *pINChI, INChI_Stereo *Stereo )
{
    int nRet = SP3_NONE;
    if ( pINChI && !pINChI->bDeleted && Stereo && 0 < Stereo->nNumberOfStereoCenters ) {
        if ( 0 != Stereo->nCompInv2Abs ) {
            if ( pINChI->nFlags & INCHI_FLAG_REL_STEREO ) {
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
                if ( 1 <  Stereo->nNumberOfStereoCenters ) {
                    nRet = SP3_REL;
                }
#else
                nRet = SP3_REL;
#endif
            } else
            if ( pINChI->nFlags & INCHI_FLAG_RAC_STEREO ) {
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
                if ( 1 <  Stereo->nNumberOfStereoCenters ) {
                    nRet = SP3_REL;
                }
#else
                nRet = SP3_RAC;
#endif
            } else {
                nRet = SP3_ABS;
            }
        } else
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
        if ( !(( pINChI->nFlags & (INCHI_FLAG_REL_STEREO|INCHI_FLAG_RAC_STEREO) ) && 1 == Stereo->nNumberOfStereoCenters) )
#endif
        {
            nRet = SP3_ONLY; /*  SP3_NONE if relative stereo and 1 stereocenter */
        }
    }
    return nRet;
}

/* char sDifSegs[DIFL_LENGTH][DIFS_LENGTH]; */

/**********************************************************************************************/
/*  sorting in descending order: return -1 if *p1 > *p2, return +1 if *p1 < *p2               */
/**********************************************************************************************/
int CompINChILayers(const INCHI_SORT *p1, const INCHI_SORT *p2, 
                    char sDifSegs[][DIFS_LENGTH], int bFixTranspChargeBug )
{
    int ret = 0, num, i, num_H1, num_H2;
    
    const INChI *i1  = NULL; /* Mobile-H layers in Mobile-H sorting order */
    const INChI *i2  = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    
    int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    
    INChI_Stereo *Stereo1, *Stereo2;
    INChI_Stereo *IsoStereo1, *IsoStereo2;
    int bRelRac[DIFL_LENGTH];
    char *psDifSegs;

    n1 = ( p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms )? TAUT_YES : TAUT_NON;

    i1  = p1->pINChI[n1];
    i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
           p2->pINChI[TAUT_NON]->nNumberOfAtoms)? p2->pINChI[TAUT_NON] : (const INChI *)NULL;

    num_H1 = num_H2 = 0;
    memset( bRelRac, DIFV_BOTH_EMPTY, sizeof(bRelRac) );
    /*=====================*/
    /*====     /f    ======*/
    /*=====================*/
    if ( i1 && !i1->bDeleted && i1->szHillFormula && i1->szHillFormula[0] ) {
        sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
        if ( i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0] ) {
            if ( !CompareHillFormulasNoH( i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2 ) &&
                  num_H1 == num_H2 ) {
                sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_EQL2PRECED;
            } else {
                sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
            }
        } else {
            sDifSegs[DIFL_F][DIFS_f_FORMULA] |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
        }
    } else {
        sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
        if ( i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0] ) {
            sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
        } else {
            sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
        }
    }
    /*=====================*/
    /*====     /c    ======*/
    /*=====================*/
    if ( i1 && !i1->bDeleted && i1->lenConnTable > 1 ) {
        sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    } else {
        sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    }
    /*=====================*/
    /*====     /h    ======*/
    /*=====================*/
    /* M: H atoms */
    if ( i1 && !i1->bDeleted ) {
        num_H1 = (i1->lenTautomer > 0 && i1->nTautomer && i1->nTautomer[0])? 1 : 0; /* number of t-groups */
        if ( !num_H1 && i1->nNum_H ) {
            for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) { /* immobile H */
                if ( i1->nNum_H[i] ) {
                    num_H1 = 1;
                    break;
                }
            }
        }
        sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= num_H1? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    } else {
        sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    }
    /* F: fixed mobile H */
    if ( i2 && !i2->bDeleted && i2->nNum_H_fixed ) {
        num_H2 = 0;
        if ( i1 && !i1->bDeleted ) {
            for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) {
                if ( i2->nNum_H_fixed[i] ) {
                    num_H2 = 1;
                    break;
                }
            }
        }
        sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= num_H2? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    } else {
        sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    }
    /* MI: exchangable isotopic H: see OutputINChI1(), num_iso_H[] */

    /*=====================*/
    /*====     /q    ======*/
    /*=====================*/
    psDifSegs = &sDifSegs[DIFL_F][DIFS_q_CHARGE];
    if ( i1 && !i1->bDeleted ) {
        if ( i1->nTotalCharge ) {
            sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
        } else {
            sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
        }
        if ( i2 && !i2->bDeleted ) {
            if ( i1->nTotalCharge ) {
                if ( i1->nTotalCharge == i2->nTotalCharge ) {
                    *psDifSegs |= DIFV_EQL2PRECED;
                } else
                if ( i2->nTotalCharge ) {
                    *psDifSegs |= DIFV_NEQ2PRECED;
                } else {
                    *psDifSegs |= DIFV_IS_EMPTY;
                }
            } else {
                if ( i2->nTotalCharge ) {
                    *psDifSegs |= DIFV_NEQ2PRECED;
                } else {
                    *psDifSegs |= DIFV_BOTH_EMPTY;
                }
            }
        } else
        if ( !i2 ) {
            if (bFixTranspChargeBug==1) 
            {
            /* bug explanation:

            component #1 is tautomeric, component #2 is not
            Mobile-H(#2) > Mobile-H(#1)
            Fixed-H(#2) = Mobile-H(#2) < Fixed-H(#1)

            Layer       first_charge   second_charge

            Mobile-H    0    (comp#1)  -1 (comp#2)
            Fixed-H     none (comp#2)  -1 (comp#1)

            v1.01 charge compared decided that charge layers are same and omitted Fixed-H /q layer

            Solution: when component permutation is detected AND fixed-H component does not exist,
            compare Mobile-H charge [0 (comp#1) in the example] to the charge of Mobile-H [-1 (comp#2)]
            of the component that has none Fixed-H charge
            */

            /* Fixed-H i2 is empty because Fixed-H struct is same as Mobile-H */
                if ( p1->ord_number != p2->ord_number && /* component order in Fixed-H is different from Mobile-H */
                     n1 == TAUT_YES && p2->pINChI[TAUT_YES] && !p2->pINChI[TAUT_YES]->bDeleted &&
                    p2->pINChI[TAUT_YES]->nNumberOfAtoms ) {
                    int i2_nTotalCharge = p2->pINChI[TAUT_YES]->nTotalCharge;

                    if ( i1->nTotalCharge ) {
                        if ( i1->nTotalCharge == i2_nTotalCharge ) {
                            *psDifSegs |= DIFV_EQL2PRECED;
                        } else
                        if ( i2_nTotalCharge ) {
                            *psDifSegs |= DIFV_NEQ2PRECED;
                        } else {
                            *psDifSegs |= DIFV_IS_EMPTY;
                        }
                    } else {
                        if ( i2_nTotalCharge ) {
                            *psDifSegs |= DIFV_NEQ2PRECED;
                        } else {
                            *psDifSegs |= DIFV_BOTH_EMPTY;
                        }
                    }
                } else {
                    *psDifSegs |= i1->nTotalCharge? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
                }
            } 
            else /* if (bFixTranspChargeBug==1) */
            {
                *psDifSegs |= i1->nTotalCharge? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
            } 
        } 

        else /* if ( !i2 ) { */
        {
            /* i2 && i2->bDeleted */
            *psDifSegs |= i1->nTotalCharge? DIFV_IS_EMPTY : DIFV_BOTH_EMPTY;
        }

    } else {
        sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
        if ( i2 && !i2->bDeleted ) {
            if ( i2->nTotalCharge ) {
                sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
            } else {
                sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
            }
        }
    }
    /*************** stereo *****************/
    if ( i1 && !i1->bDeleted ) {
        Stereo1    = i1->Stereo;
        IsoStereo1 = i1->StereoIsotopic;
    } else {
        Stereo1    = NULL;
        IsoStereo1 = NULL;
    }
    if ( i2 && !i2->bDeleted ) {
        Stereo2    = i2->Stereo;
        IsoStereo2 = i2->StereoIsotopic;
    } else {
        Stereo2    = NULL;
        IsoStereo2 = NULL;
    }
    /*=====================*/
    /*====     /b    ======*/
    /*=====================*/
    /* M double bond stereo */
    psDifSegs = &sDifSegs[DIFL_M][DIFS_b_SBONDS];
    if ( Stereo1 && Stereo1->nNumberOfStereoBonds ) {
        *psDifSegs |= DIFV_NEQ2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* F double bond stereo */
    psDifSegs = &sDifSegs[DIFL_F][DIFS_b_SBONDS];
    if ( Stereo2 && Stereo2->nNumberOfStereoBonds ) {
        if ( Stereo1 && Stereo1->nNumberOfStereoBonds ) {
            if ( Eql_INChI_Stereo( Stereo1, EQL_SP2, Stereo2, EQL_SP2, 0 ) ) {
                *psDifSegs |= DIFV_EQL2PRECED;
            } else {
                *psDifSegs |= DIFV_NEQ2PRECED;
            }
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else {
        if ( Stereo1 && Stereo1->nNumberOfStereoBonds ) {
            *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_BOTH_EMPTY;
        }
    }
    /* MI double bond stereo */
    psDifSegs = &sDifSegs[DIFL_MI][DIFS_b_SBONDS];
    if ( IsoStereo1 && IsoStereo1->nNumberOfStereoBonds ) {
        if ( Eql_INChI_Stereo( IsoStereo1, EQL_SP2, Stereo1, EQL_SP2, 0 ) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else {
        if ( Stereo1 && Stereo1->nNumberOfStereoBonds ) {
            *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
        } else {
            *psDifSegs |= DIFV_BOTH_EMPTY;
        }
    }
    /* FI double bond stereo */
    psDifSegs = &sDifSegs[DIFL_FI][DIFS_b_SBONDS];
    if ( IsoStereo2 && IsoStereo2->nNumberOfStereoBonds ) {
        if ( Eql_INChI_Stereo( IsoStereo2, EQL_SP2, Stereo2, EQL_SP2, 0 ) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            if ( !(Stereo1 && Stereo1->nNumberOfStereoBonds) &&
                 !(Stereo2 && Stereo2->nNumberOfStereoBonds) &&
                 Eql_INChI_Stereo( IsoStereo2, EQL_SP2, IsoStereo1, EQL_SP2, 0 ) ) {
                *psDifSegs |= DIFV_FI_EQ_MI;
            } else {
                *psDifSegs |= DIFV_NEQ2PRECED;
            }
        }
    } else {
        /* the solution table for FI stereo,
           in case of FI stereo is empty
           E = segment is empty, NE = not empty
           +==============================+
           | M   | MI  | F   |  result    |
           +=====+=====+=====+============+
           | E   | E   | E   | both empty |
           +-----+-----+-----+------------+
           | NE  | E   | E   | both empty |
           +-----+-----+-----+------------+
           | E   | NE  | E   | is empty   |
           +-----+-----+-----+------------+
           | NE  | NE  | E   | both empty |
           +-----+-----+-----+------------+
           | E   | E   | NE  | is empty   |
           +-----+-----+-----+------------+
           | NE  | E   | NE  | is empty   |
           +-----+-----+-----+------------+
           | E   | NE  | NE  | is empty   |
           +-----+-----+-----+------------+
           | NE  | NE  | ME  | is empty   |
           +==============================+
        */
        if ( Stereo2    && Stereo2->nNumberOfStereoBonds ) {
            *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
        } else
        if ( IsoStereo1 && IsoStereo1->nNumberOfStereoBonds &&
             !(Stereo1 && Stereo1->nNumberOfStereoBonds)
            ) {
            *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_BOTH_EMPTY;
        }
    }
    /*==================================*/
    /*====     /t, /m, /s for M   ======*/
    /*==================================*/
    /* M sp3 stereo */


    bRelRac[DIFL_M ] = GetSp3RelRacAbs( i1, Stereo1 );       /* Mobile-H */
    bRelRac[DIFL_MI] = GetSp3RelRacAbs( i1, IsoStereo1 );
    bRelRac[DIFL_F ] = GetSp3RelRacAbs( i2, Stereo2 );       /* Fixed-H */
    bRelRac[DIFL_FI] = GetSp3RelRacAbs( i2, IsoStereo2 );
    if ( SP3_NONE != bRelRac[DIFL_M] ) {
        sDifSegs[DIFL_M][DIFS_t_SATOMS] |= (bRelRac[DIFL_M] & SP3_ANY)?  DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
        sDifSegs[DIFL_M][DIFS_m_SP3INV] |= (bRelRac[DIFL_M] & SP3_ABS)?  DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
        sDifSegs[DIFL_M][DIFS_s_STYPE]  |= (bRelRac[DIFL_M] & SP3_TYPE)? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    } else {
        sDifSegs[DIFL_M][DIFS_t_SATOMS] |= DIFV_BOTH_EMPTY;
        sDifSegs[DIFL_M][DIFS_m_SP3INV] |= DIFV_BOTH_EMPTY;
        sDifSegs[DIFL_M][DIFS_s_STYPE]  |= DIFV_BOTH_EMPTY;
    }
    /*=====================*/
    /*====     /t    ======*/
    /*=====================*/
    /* F sp3 stereo */
    psDifSegs = &sDifSegs[DIFL_F][DIFS_t_SATOMS];
    if ( SP3_ANY & bRelRac[DIFL_F] ) {
        if ( Eql_INChI_Stereo( Stereo2, EQL_SP3, Stereo1, EQL_SP3, 0 ) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else
    if ( SP3_ANY & bRelRac[DIFL_M] ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* MI sp3 stereo */
    psDifSegs = &sDifSegs[DIFL_MI][DIFS_t_SATOMS];
    if ( SP3_ANY & bRelRac[DIFL_MI] ) {
        if ( Eql_INChI_Stereo( IsoStereo1, EQL_SP3, Stereo1, EQL_SP3, 0 ) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else
    if ( SP3_ANY & bRelRac[DIFL_M] ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* FI sp3 stereo */
    psDifSegs = &sDifSegs[DIFL_FI][DIFS_t_SATOMS];
    if ( SP3_ANY & bRelRac[DIFL_FI] ) {
        if ( Eql_INChI_Stereo( IsoStereo2, EQL_SP3, Stereo2, EQL_SP3, 0 ) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else
        if ( !(SP3_ANY & bRelRac[DIFL_M]) &&
             !(SP3_ANY & bRelRac[DIFL_F]) &&
             Eql_INChI_Stereo( IsoStereo2, EQL_SP3, IsoStereo1, EQL_SP3, 0 ) ) {
            *psDifSegs |= DIFV_FI_EQ_MI;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else /* similar to /b */
    if ( (SP3_ANY & bRelRac[DIFL_F]) ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else
    if ( (SP3_ANY & bRelRac[DIFL_MI]) && !(SP3_ANY & bRelRac[DIFL_M]) ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /*=====================*/
    /*====     /m    ======*/
    /*=====================*/
    /* F sp3 abs stereo inversion */
    psDifSegs = &sDifSegs[DIFL_F][DIFS_m_SP3INV];
    if ( bRelRac[DIFL_F] & SP3_ABS ) {
        /* the order of || operands below is critically important: || is not a commutative operation */
        if ( !(bRelRac[DIFL_M] & SP3_ABS) || Stereo2->nCompInv2Abs != Stereo1->nCompInv2Abs ) {
            *psDifSegs |= DIFV_NEQ2PRECED;
        } else {
            *psDifSegs |= DIFV_EQL2PRECED;
        }
    } else
    if ( bRelRac[DIFL_M] & SP3_ABS ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* MI sp3 abs stereo inversion */
    psDifSegs = &sDifSegs[DIFL_MI][DIFS_m_SP3INV];
    if ( SP3_ABS & bRelRac[DIFL_MI] ) {
        if ( (SP3_ABS & bRelRac[DIFL_M]) && IsoStereo1->nCompInv2Abs == Stereo1->nCompInv2Abs ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else
    if ( SP3_ABS & bRelRac[DIFL_M] ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* FI sp3 abs stereo inversion */
    psDifSegs = &sDifSegs[DIFL_FI][DIFS_m_SP3INV];
    if ( SP3_ABS & bRelRac[DIFL_FI] ) {
        if ( (SP3_ABS & bRelRac[DIFL_F]) && IsoStereo2->nCompInv2Abs == Stereo2->nCompInv2Abs ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else
        if ( !(SP3_ABS & bRelRac[DIFL_M]) &&
             !(SP3_ABS & bRelRac[DIFL_F]) &&
             (SP3_ABS & bRelRac[DIFL_MI]) && /* make sure IsoStereo1 != NULL */
             IsoStereo2->nCompInv2Abs == IsoStereo1->nCompInv2Abs ) {
            *psDifSegs |= DIFV_FI_EQ_MI;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else /* similar to /b */
    /* the order of || operands below is critically important: || is no a commutative operation */
    if ( (SP3_ABS & bRelRac[DIFL_F]) ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else
    if ( (SP3_ABS & bRelRac[DIFL_MI]) && !(SP3_ABS & bRelRac[DIFL_M]) ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /*=====================*/
    /*====     /s    ======*/
    /*=====================*/
    /* F sp3 stereo type */
    psDifSegs = &sDifSegs[DIFL_F][DIFS_s_STYPE];
    if ( bRelRac[DIFL_F] & SP3_TYPE ) {
        if ( (bRelRac[DIFL_F] & SP3_TYPE) == (bRelRac[DIFL_M] & SP3_TYPE) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else
    if ( bRelRac[DIFL_M] & SP3_TYPE ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* MI sp3 stereo type */
    psDifSegs = &sDifSegs[DIFL_MI][DIFS_s_STYPE];
    if ( SP3_TYPE & bRelRac[DIFL_MI] ) {
        if ( (SP3_TYPE & bRelRac[DIFL_MI]) == (SP3_TYPE & bRelRac[DIFL_M]) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else
    if ( SP3_TYPE & bRelRac[DIFL_M] ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* FI sp3 stereo type */
    psDifSegs = &sDifSegs[DIFL_FI][DIFS_s_STYPE];
    if ( SP3_TYPE & bRelRac[DIFL_FI] ) {
        if ( (SP3_TYPE & bRelRac[DIFL_FI]) == (SP3_TYPE & bRelRac[DIFL_F]) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else
        if ( !(SP3_TYPE & bRelRac[DIFL_M]) &&
             !(SP3_TYPE & bRelRac[DIFL_F]) &&
             (SP3_TYPE & bRelRac[DIFL_MI]) ) {
            *psDifSegs |= DIFV_FI_EQ_MI;
        } else {
            *psDifSegs |= DIFV_NEQ2PRECED;
        }
    } else /* similar to /b */
    /* the order of || operands below is critically important: || is not a commutative operation */
    if ( (SP3_TYPE & bRelRac[DIFL_F]) ) {
        *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    } else
    if ( (SP3_TYPE & bRelRac[DIFL_MI]) && !(SP3_TYPE & bRelRac[DIFL_M]) ) {
        *psDifSegs |= i2? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /*=====================*/
    /*====     /o    ======*/
    /*=====================*/
    if ( p1 && p2 && p1->ord_number != p2->ord_number ) {
        sDifSegs[DIFL_F][DIFS_o_TRANSP] |= DIFV_NEQ2PRECED;
    }
    /*=====================*/
    /*====     /i    ======*/
    /*=====================*/
    /* M isotopic atoms */
    psDifSegs = &sDifSegs[DIFL_MI][DIFS_i_IATOMS];
    if ( i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups) ) {
        *psDifSegs |= DIFV_NEQ2PRECED;
    } else {
        *psDifSegs |= DIFV_BOTH_EMPTY;
    }
    /* F isotopic atoms */
    psDifSegs = &sDifSegs[DIFL_FI][DIFS_i_IATOMS];
    if ( i2 && !i2->bDeleted ) {
        if ( i2->nNumberOfIsotopicAtoms || i2->nNumberOfIsotopicTGroups ) {
            if ( !i1 || i1->bDeleted ||
                 i2->nNumberOfIsotopicAtoms   != i1->nNumberOfIsotopicAtoms ||
                 i2->nNumberOfIsotopicTGroups != i1->nNumberOfIsotopicTGroups ) {
                *psDifSegs |= DIFV_NEQ2PRECED;
            } else {
                int diff;
                num  = i1->nNumberOfIsotopicAtoms;
                diff = 0;
                for ( i = 0; i < num; i ++ ) {
                    /* compare isotopic atoms */
                    if ( diff = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber )
                        break;
                    if ( diff = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference )
                        break;
                    /* compare isotopic H */
                    if ( diff = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T )
                        break;
                    if ( diff = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D )
                        break;
                    if ( diff = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H )
                        break;
                }
                if ( !diff ) {
                    num = i1->nNumberOfIsotopicTGroups;
                    for ( i = 0; i < num; i ++ ) {
                        if ( diff = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber )
                            break;
                        if ( diff = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T )
                            break;
                        if ( diff = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D )
                            return diff;
                        if ( diff = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H )
                            break;
                    }
                }
                *psDifSegs |= diff? DIFV_NEQ2PRECED : DIFV_FI_EQ_MI;

            }
        } else
        if ( i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups) ) {
            *psDifSegs |= DIFV_IS_EMPTY;
        }
    } else
    if ( !i2 ) {
        if ( i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups) ) {
            *psDifSegs |= DIFV_EQL2PRECED;
        } else {
            *psDifSegs |= DIFV_BOTH_EMPTY;
        }
    }

    return ret;
}
/**********************************************************************************************/
int INChI_SegmentAction( char cDifSegs )
{
    if ( !(cDifSegs & DIFV_OUTPUT_OMIT_F) ) {
        return INCHI_SEGM_OMIT;
    }
    if ( (cDifSegs & DIFV_OUTPUT_EMPTY_T) && !(cDifSegs & DIFV_OUTPUT_EMPTY_F) ) {
        return INCHI_SEGM_EMPTY;
    }
    if ( (cDifSegs & DIFV_OUTPUT_FILL_T) ) {
        return INCHI_SEGM_FILL;
    }
    return INCHI_SEGM_OMIT; /* the control flow shoul never reach this point */
}
/**********************************************************************************************/
int MarkUnusedAndEmptyLayers( char sDifSegs[][DIFS_LENGTH] )
{
    /****************************************************
     1. If all elements of a layer are DIFV_IS_EMPTY and/or DIFV_BOTH_EMPTY
        and/or DIFV_EQL2PRECED and/or DIFV_FI_EQ_MI
        and there is NO succeeding non-empty layer then mark the 1st element
        of the layer DIFV_BOTH_EMPTY; this layerr will be omitted.
     
     2. If all elements of a layer are DIFV_IS_EMPTY and/or DIFV_BOTH_EMPTY
        and/or DIFV_EQL2PRECED and/or DIFV_FI_EQ_MI
        and there IS a succeeding non-empty layer then mark the 1st element
        of the layer DIFV_IS_EMPTY and all other elements DIFV_BOTH_EMPTY;
        only the first empty segment of this layerr will be output.

     3. If NOT all elements of a layer are DIFV_IS_EMPTY and/or DIFV_BOTH_EMPTY
        and/or DIFV_EQL2PRECED and/or DIFV_FI_EQ_MI
        and the 1st element of the layer is DIFV_BOTH_EMPTY then mark it
        DIFV_IS_EMPTY; it will be output as empty (except M layer).
     */

    int i, nLayer, sBits, nFirstSegm;
#define nFirstFmlSegm   DIFS_f_FORMULA
#define nFirstIsoSegm   DIFS_i_IATOMS
     /* FI */
    nLayer = DIFL_FI;
    nFirstSegm = nFirstIsoSegm;
    sBits  = 0;
    for ( i = 0; i < DIFS_idf_LENGTH; i ++ ) {
        sBits |= sDifSegs[nLayer][i];
    }
    if ( !(sBits & DIFV_OUTPUT_OMIT_F) ) {
        /* Omit the FI layer */
        memset( sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH);
    } else
    if ( sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
         !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F) ) {
        sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    }

    /* MI */
    nLayer = DIFL_MI;
    nFirstSegm = nFirstIsoSegm;
    sBits  = 0;
    for ( i = 0; i < DIFS_idf_LENGTH; i ++ ) {
        sBits |= sDifSegs[nLayer][i];
    }
    if ( !(sBits & DIFV_OUTPUT_OMIT_F) ) {
        /* Omit the MI layer */
        memset( sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH);
    } else
    if ( sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
         !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F) ) {
        sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    }

    /* F */
    nLayer = DIFL_F;
    nFirstSegm = nFirstFmlSegm;
    sBits  = 0;
    for ( i = 0; i < DIFS_idf_LENGTH; i ++ ) {
        sBits |= sDifSegs[nLayer][i];
    }
    if ( !(sBits & DIFV_OUTPUT_OMIT_F) &&
         sDifSegs[DIFL_FI][nFirstIsoSegm] == DIFV_BOTH_EMPTY ) {
        /* Omit the F layer: no non-iotopic and no isotopic segments */
        memset( sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH);
    } else
    /* do not omit fixed-H layer */
    if ( sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
         !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F) ) {
        sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    }
     
    /* M -- leave as it is */
    return 0;
#undef nFirstFmlSegm
#undef nFirstIsoSegm
}
/*********************************************************************************************/
int CompareInchiStereo( INChI_Stereo *Stereo1, INCHI_MODE nFlags1, INChI_Stereo *Stereo2, INCHI_MODE nFlags2 )
{
    int i, num, ret;
    if ( Stereo2 && Stereo1 ) {
        /*  compare stereogenic bonds */
        num = inchi_min( Stereo2->nNumberOfStereoBonds, Stereo1->nNumberOfStereoBonds );
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)Stereo2->nBondAtom1[i] - (int)Stereo1->nBondAtom1[i] )
                return ret;
            if ( ret = (int)Stereo2->nBondAtom2[i] - (int)Stereo1->nBondAtom2[i] )
                return ret;
            if ( ret = (int)Stereo2->b_parity[i] - (int)Stereo1->b_parity[i] )
                return ret;
        }
        if ( ret = (int)Stereo2->nNumberOfStereoBonds - (int)Stereo1->nNumberOfStereoBonds )
            return ret;
        /*  compare stereogenic atoms */
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
        if ( ((nFlags1 | nFlags2) & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO) ) &&
             1 == Stereo2->nNumberOfStereoCenters &&
             1 == Stereo1->nNumberOfStereoCenters ) {
            ; /*  do not compare single stereocenters in case of relative stereo */
        } else
#endif
        {
            num = inchi_min( Stereo2->nNumberOfStereoCenters, Stereo1->nNumberOfStereoCenters );
            for ( i = 0; i < num; i ++ ) {
                if ( ret = (int)Stereo2->nNumber[i] - (int)Stereo1->nNumber[i] )
                    return ret;
                if ( ret = (int)Stereo2->t_parity[i] - (int)Stereo1->t_parity[i] )
                    return ret;
            }
            if ( ret = (int)Stereo2->nNumberOfStereoCenters - (int)Stereo1->nNumberOfStereoCenters )
                return ret;
            /*  compare stereo-abs-is-inverted flags  for non-relative, non-racemic */
            if ( !((nFlags1 | nFlags2) & (INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO)) ) {
                if ( ret = (Stereo2->nCompInv2Abs < 0) - (Stereo1->nCompInv2Abs < 0) ) {
                    return ret;
                }
            }
        }
    } else
    if ( Stereo2 && ( Stereo2->nNumberOfStereoBonds > 0 ||
                      Stereo2->nNumberOfStereoCenters > 0
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
                      && /*  do not compare single stereocenters in case of relative stereo */
                      !((nFlags2 & (INCHI_FLAG_REL_STEREO|INCHI_FLAG_RAC_STEREO)) &&
                       1 == Stereo2->nNumberOfStereoCenters
                       )
#endif
       ) ) {
        return 1;
    }else
    if ( Stereo1 && ( Stereo1->nNumberOfStereoBonds > 0 || 
                      Stereo1->nNumberOfStereoCenters > 0
#if( REL_RAC_STEREO_IGN_1_SC == 1 )
                      && /*  do not compare single stereocenters in case of relative stereo */
                      !((nFlags1 & (INCHI_FLAG_REL_STEREO|INCHI_FLAG_RAC_STEREO)) &&
                       1 == Stereo1->nNumberOfStereoCenters
                       )
#endif
       ) ) {
        return -1;
    }
    return 0;
}
/**********************************************************************************************/
/*  sorting in descending order: return -1 if *p1 > *p2, return +1 if *p1 < *p2               */
/**********************************************************************************************/
int CompINChI2(const INCHI_SORT *p1, const INCHI_SORT *p2, int bTaut, int bCompareIsotopic)
{
    int ret, num, i, num_H1, num_H2;
    
    const INChI *i1  = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    const INChI *i2  = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    
    int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    int   n2;               /* TAUT_YES if tautomeric i2 exists, otherwise TAUT_NON */
    
    const INChI *i1n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */
    const INChI *i2n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */

    /*const INChI *i1t = NULL;*/ /* temp for i1 if both tautomeric AND non-tautomeric exist */
    /*const INChI *i2t = NULL;*/ /* temp for i2 if both tautomeric AND non-tautomeric exist */


    /* INChI_Stereo *Stereo1, *Stereo2; */

    n1 = ( p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms )? TAUT_YES : TAUT_NON;
    n2 = ( p2->pINChI[TAUT_YES] && p2->pINChI[TAUT_YES]->nNumberOfAtoms )? TAUT_YES : TAUT_NON;

    i1  = p1->pINChI[n1];
    i1n = (n1 == TAUT_YES && p1->pINChI[TAUT_NON] &&
           p1->pINChI[TAUT_NON]->nNumberOfAtoms)? p1->pINChI[TAUT_NON] : (const INChI *)NULL;

    i2  = p2->pINChI[n2];
    i2n = (n2 == TAUT_YES && p2->pINChI[TAUT_NON] &&
          p2->pINChI[TAUT_NON]->nNumberOfAtoms)? p2->pINChI[TAUT_NON] : (const INChI *)NULL;

    /* non-deleted-non-empty < deleted < empty */
    if ( i1 && !i2 )
        return -1;   /* non-empty is the smallest (first) */
    if ( !i1 && i2 )
        return 1;
    if ( !i1 && !i2 )
        return 0;
    if ( i1->bDeleted && !i2->bDeleted )
        return 1;    /* deleted is the largest (last) among non-empty */
    if ( !i1->bDeleted && i2->bDeleted )
        return -1;

    num_H1 = num_H2 = 0;
    
    /* do not compare terminal H */
    if ( ret = CompareHillFormulasNoH( i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2 ) ) {
        return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    }

    /*********************************************************
            compare non-isotopic non-tautomeric part
     *********************************************************/

    /* compare number of atoms (excluding terminal H) */
    if ( ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms )
        return ret; /*  more atoms first */
    
    /*  compare elements  (excluding terminal H) */
    num = i1->nNumberOfAtoms;
    for ( i = 0; i < num; i ++ ) { /* should always be equal if Hill formulas are same */
        if ( ret = (int)i2->nAtom[i] - (int)i1->nAtom[i] )
            return ret; /* greater periodic number first */
    }
    /**********************************************************
        compare connection tables
    ***********************************************************/
    if ( ret = i2->lenConnTable - i1->lenConnTable )
        return ret; /* longer connection table first */
    num = i2->lenConnTable; 
    for ( i = 0; i < num; i ++ ) {
        if ( ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i] )
            return ret; /* greater connection table first */
    }
    /*********************************************************
      compare compare total number of H (inverse: H3 < H2 )
    **********************************************************/
    if ( ret = num_H2 - num_H1 )
        return ret;
    /*********************************************************
      compare non-tautomeric num_H: N < NH3 < NH2 < NH
    **********************************************************/
    num = i1->nNumberOfAtoms;
    for ( i = 0; i < num; i ++ ) {
        if ( i2->nNum_H[i] != i1->nNum_H[i] ) {
            return !i2->nNum_H[i]?  1 :  /* no H first */
                   !i1->nNum_H[i]? -1 :
                   (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
        }
    }
    /*********************************************************
         compare non-isotopic tautomeric part
     *********************************************************/
    if ( ret = CompareTautNonIsoPartOfINChI( i1, i2) ) {
        return ret;
    }
    /*
    if ( ret = i2->lenTautomer - i1->lenTautomer )
        return ret;
    num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    for ( i = 0; i < num; i ++ ) {
        if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
            return ret;
    }
    */
    /*********************************************************
     *                                                       *
     *  at this point both components are either tautomeric  *
     *  or non-tautomeric                                    *
     *                                                       *
     *********************************************************/

    /*********************************************************
       non-tautomeric "fixed H" specific
     *********************************************************/
    if ( TAUT_NON == bTaut && (i1n && i1n->nNum_H_fixed || i2n && i2n->nNum_H_fixed) ) {
        /* first, compare non-tautomeric chem. formulas -- they may be different */
        const char *f1 = (i1n /*&& i1n->nNum_H_fixed*/)? i1n->szHillFormula : i1->szHillFormula;
        const char *f2 = (i2n /*&& i2n->nNum_H_fixed*/)? i2n->szHillFormula : i2->szHillFormula;
        if ( f1 && f2 &&(ret = CompareHillFormulas( f1, f2 ))) {
            return ret;
        }
        /* secondly, compare fixed-H distribution */
        if ( i1n && i1n->nNum_H_fixed && i2n && i2n->nNum_H_fixed ) {
            num = inchi_min( i1n->nNumberOfAtoms, i2n->nNumberOfAtoms);
            for ( i = 0; i < num; i ++ ) {
                if ( i2n->nNum_H_fixed[i] != i1n->nNum_H_fixed[i] ) {
                    return !i2n->nNum_H_fixed[i]?  1 : /* no fixed H first */
                           !i1n->nNum_H_fixed[i]? -1 :
                           (int)i2n->nNum_H_fixed[i] - (int)i1n->nNum_H_fixed[i];
                }
            }
            if ( ret = (int)i2n->nNumberOfAtoms - (int)i1n->nNumberOfAtoms ) {
                return ret; /* should not happen <BRKPT> */
            }
        } else
        if ( i1n && i1n->nNum_H_fixed ) {
            num = i1n->nNumberOfAtoms;
            for ( i = 0; i < num; i ++ ) {  /* added 2004-05-04 */
                if ( i1n->nNum_H_fixed[i] ) {
                    return -1; /* i1n->nNum_H_fixed[i] > 0? -1:1;*/
                }
            }
            /* p1 is tautomeric, p2 is not tautomeric; this must have been detected earlier */
            /*return -1;*/ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
        } else {
            num = i2n->nNumberOfAtoms;
            for ( i = 0; i < num; i ++ ) {  /* added 2004-05-04 */
                if ( i2n->nNum_H_fixed[i] ) {
                    return 1; /* i2n->nNum_H_fixed[i] > 0? 1:-1;*/
                }
            }
            /* p2 is tautomeric, p1 is not tautomeric; this must have been detected earlier */
            /*return 1; */ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
        }
    }
    
    /*************************************************************************
              if requested non-tautomeric comparison then
              prepare to compare non-taut non-isotopic stereo, etc. 
     *************************************************************************/
    if ( TAUT_NON == bTaut ) {
        if ( i1n ) {
            /*i1t = i1;*/
            i1  = i1n;
        }
        if ( i2n ) {
            /*i2t = i2;*/
            i2  = i2n;
        }
    }

    /*********************************************************
        compare non-isotopic stereo
     *********************************************************/
    ret = CompareInchiStereo( i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags );
    if ( ret ) {
        return ret;
    }
    /*******************************************************
        do not switch back to tautomeric i1, i2
     *******************************************************/
    /* -- how to switch back --
    if ( i1t ) {
        i1  = i1t;
        i1t = NULL;
    }
    if ( i2t ) {
        i2  = i2t;
        i2t = NULL;
    }
    */
    /******************************************************
         compare isotopic non-tautomeric part
     ******************************************************/
    if ( bCompareIsotopic ) {
        if ( ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms )
            return ret;
        num = i1->nNumberOfIsotopicAtoms;
        /*  compare isotopic atoms */
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference )
                return ret;
        }
        /* compare isotopic H */
        /* if tautomeric comparison mode then here are compared only non-tautomeric H */
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D )
                return ret;
            if ( ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H )
                return ret;
        }
        /*****************************************************
             compare isotopic tautomeric part
         *****************************************************/
        if ( ret = i2->nNumberOfIsotopicTGroups - i1->nNumberOfIsotopicTGroups )
            return ret;
        num = i1->nNumberOfIsotopicTGroups;
        for ( i = 0; i < num; i ++ ) {
            if ( ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D )
                return ret;
            if ( ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H )
                return ret;
        }
    
        /****************************************************
            compare isotopic stereo
         ****************************************************/
        ret = CompareInchiStereo( i1->StereoIsotopic, i1->nFlags, i2->StereoIsotopic, i2->nFlags );
        if ( ret ) {
            return ret;
        }
    }

    /**********************************************************
        compare charges: non-charged first, then in order of
        ascending charges (negative first) 
    ***********************************************************/
    if ( i2->nTotalCharge && i1->nTotalCharge ) {
        /*  both are charged; smaller charges first */
        ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
        return ret;
    }
    if ( ret = (i1->nTotalCharge? 1:0) - (i2->nTotalCharge? 1:0) ) {
        /*  only one is charged; uncharged first */
        return ret;
    }
    /* stable sort */
    /*ret = p1->ord_number - p2->ord_number;*/

    return ret;
}
/***********************************************************************/
int CompINChINonTaut2(const void *p1, const void *p2)
{
    int ret;
    ret = CompINChI2( (const INCHI_SORT *)p1, (const INCHI_SORT *)p2, TAUT_NON, 1 );
#if( CANON_FIXH_TRANS == 1 )
    if ( !ret ) {
        /* to obtain canonical transposition 2004-05-10 */
        ret = CompINChI2( (const INCHI_SORT *)p1, (const INCHI_SORT *)p2, TAUT_YES, 1 );
    }
#endif
    if ( !ret ) {
        /* stable sort */
        ret = ((const INCHI_SORT *)p1)->ord_number - ((const INCHI_SORT *)p2)->ord_number;
    }
    return ret;
}
/***********************************************************************/
int CompINChITaut2(const void *p1, const void *p2)
{
    int ret;
    ret = CompINChI2( (const INCHI_SORT *)p1, (const INCHI_SORT *)p2, TAUT_YES, 1 );
#if( CANON_FIXH_TRANS == 1 )
    if ( !ret ) {
        /* to obtain canonical transposition 2004-05-10 */
        ret = CompINChI2( (const INCHI_SORT *)p1, (const INCHI_SORT *)p2, TAUT_NON, 1 );
    }
#endif
    if ( !ret ) {
        /* stable sort */
        ret = ((const INCHI_SORT *)p1)->ord_number - ((const INCHI_SORT *)p2)->ord_number;
    }
    return ret;
}
/**********************************************************************************************/
/*  strrev from K&R is not in ANSI-compatible C library */
void mystrrev( char *p )
{
    char c, *q = p;
    while( *q++ )
        ;
    q -= 2; /*  pointer to the last character */
    while ( p < q ) {
        c    = *q;  /*  swap */
        *q-- = *p;
        *p++ = c;
    }
}
/*****************************************************************************************/
/*                Find DFS order for CT(canon. numbers and Hs) output                    */
/*****************************************************************************************/

static AT_NUMB   *gDfs4CT_nDfsNumber;
static AT_NUMB   *gDfs4CT_nNumDescendants;
static int        gDfs4CT_nCurrentAtom;

/**********************************************************************************************/
static int CompareDfsDescendants4CT( const void *a1, const void *a2 )
{
    int neigh1 = (int)*(const AT_RANK*)a1;
    int neigh2 = (int)*(const AT_RANK*)a2;
    if ( neigh1 > MAX_ATOMS ) {
        if ( neigh2 > MAX_ATOMS ) {
            return 0;
        }
        return 1;
    } else
    if ( neigh2 > MAX_ATOMS ) {
        return -1;
    } else {
        AT_RANK nCurDfsNumber = gDfs4CT_nDfsNumber[gDfs4CT_nCurrentAtom];
        int nDesc1 = nCurDfsNumber > gDfs4CT_nDfsNumber[neigh1]?
                         0 : (int)gDfs4CT_nNumDescendants[neigh1];
        int nDesc2 = nCurDfsNumber > gDfs4CT_nDfsNumber[neigh2]?
                         0 : (int)gDfs4CT_nNumDescendants[neigh2];
        int ret;
        if ( ret = nDesc1 - nDesc2 ) {
            return ret;
        }
        return  (int)neigh1 - (int)neigh2; /*  canon. numbers difference */
    }
}
/**********************************************************************************************/
/*  sp_ATOM *at, AT_RANK *nRank, int num_atoms */
AT_NUMB *GetDfsOrder4CT( AT_NUMB *LinearCT, int nLenCT, S_CHAR *nNum_H, int num_atoms, int nCtMode )
{
    AT_NUMB    *nStackAtom = NULL;
    int         nTopStackAtom=-1;
    AT_NUMB    *nNumDescendants = NULL; /*  number of descendants incl. closures and the atom itself */
    AT_NUMB    *nDfsNumber = NULL;
    S_CHAR     *cNeighNumb = NULL;
    NEIGH_LIST *nl = NULL;
    AT_NUMB     nDfs;
    int         i, j, u, k, start, num_rings, nTotOutputStringLen;
    AT_NUMB    *nOutputString = NULL, cDelim;
    int         bCtPredecessors = (nCtMode & CT_MODE_PREDECESSORS);

    /* int nNumStartChildren; */


    /*  allocate arrays */
    nStackAtom      = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nStackAtom[0]));
    nNumDescendants = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nNumDescendants[0]));
    nDfsNumber      = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nDfsNumber[0]));
    cNeighNumb      = (S_CHAR  *)inchi_malloc(num_atoms*sizeof(cNeighNumb[0]));
    nl              = CreateNeighListFromLinearCT( LinearCT, nLenCT, num_atoms );
    /*  check allocation */
    if ( !nStackAtom || !nNumDescendants || !nDfsNumber || !cNeighNumb || !nl ) {
        /* ret = CT_OUT_OF_RAM; */ /*  program error */ /*   <BRKPT> */
        goto exit_function;
    }
    if ( bCtPredecessors ) {
        start = 0;
    } else {
        /*  find DFS start vertex (atom) */
        for ( i = 1, start = 0; i < num_atoms; i ++ ) {
            if ( nl[i][0] <  nl[start][0] ) { /*  index = nRank-1 */
                start = i;
            }
        }
    }
    /*
      vertex information:
        1. Number of (forward edges) + (back edges, first visit -- ring closures): nl[i][0]
        2. Number of vertices traversed from this vertex, including the vertex:    nNumDescendants[i]
        3. Each edge information:
           a. forward edge (0) or back edge (1) indicator: nDfsNumber[i] > nDfsNumber[neigh]
           b. neighbor at another end of the edge neigh = nl[i][k+1], k < i

        Total per edge: 2 + 2*(number of edges)
    */

    /* DFS initiation */
    u               = start; /* start atom */
    nDfs            = 0;
    nTopStackAtom   =-1;
    memset( nDfsNumber,      0, num_atoms*sizeof(nDfsNumber[0]));
    memset( nNumDescendants, 0, num_atoms*sizeof(nNumDescendants[0]));
    memset( cNeighNumb,      0, num_atoms*sizeof(cNeighNumb[0]));
    /*  push the start atom on the stack */
    nDfsNumber[u] = ++nDfs;
    if ( bCtPredecessors ) {
        nNumDescendants[u] = 0; /* atom #1 has no predecessor */
    } else {
        nNumDescendants[u] = 1; /* count itself as a descendant */
    }
    nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    /* nNumStartChildren = 0; */
    num_rings = 0;

    /* DFS */
    
    do {
        /* advance */
        while ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i]+1,  (int)nl[i][0] >= j )
        /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
        /* replaced due to missing sequence point; undefined behavior, pointed by Geoffrey Hutchison */
        {
            cNeighNumb[i] ++;
            u = (int)nl[i][j]; /*  jth neighbor of the vertex i */
            if ( !nDfsNumber[u] ) {
                /* tree edge, 1st visit -- advance */
                /* put unexplored vertex u on the stack for further examination */
                nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
                nDfsNumber[u] = ++nDfs;
                if ( bCtPredecessors ) {
                    nNumDescendants[u] = i+1; /* predecessor's rank */
                } else {
                    nNumDescendants[u] ++; /* count atom u as its descendant */
                }
            } else
            if ( nTopStackAtom && u != (int)nStackAtom[nTopStackAtom-1] &&
                 /* back edge: u is not a predecessor of i */
                 nDfsNumber[u] < nDfsNumber[i] ) {
                /* Back edge, 1st visit: u is an ancestor of i (ring closure) */
                if ( !bCtPredecessors ) {
                    nNumDescendants[i] ++; /* count closures as descendants */
                }
                num_rings ++;          /* count ring closures */
            } else {
                nl[i][j] = MAX_ATOMS+1; /* back edge, 2nd visit: mark as deleted */
            }
        }
        cNeighNumb[i] = 0; /* all neighbors of the ith atom have been
                              traversed; resore the neighbor counter */
        /* back up */
        if ( !bCtPredecessors && nTopStackAtom /* that is, i != start */) {
            u = (int)nStackAtom[nTopStackAtom-1]; /* predecessor of i */
            nNumDescendants[u] += nNumDescendants[i]; /* add descendants */
        }
    } while ( --nTopStackAtom >= 0 );

    /* Sort the neighbors in ascending order so that:
       primary key   = number of descendants in the DFS tree; closure neighbor is 0
       secondary key = canonical number (here vertex number = canonical number - 1)
     */

    /* set static globals for the sorting: */
    gDfs4CT_nDfsNumber      = nDfsNumber;     
    gDfs4CT_nNumDescendants = nNumDescendants;
    gDfs4CT_nCurrentAtom    = -1;   

    /* sorting; deleted will be the last neighbors */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( nl[i][0] > 1 ) {
            gDfs4CT_nCurrentAtom = i;
            insertions_sort( &nl[i][1], nl[i][0], sizeof(nl[i][1]), CompareDfsDescendants4CT );
        }
        /* reduce number of neighbors to exclude deleted */
        for ( k = 0; k < nl[i][0] && nl[i][k+1] <= MAX_ATOMS; k ++ )
            ;
        nl[i][0] = k;
    }

    nTotOutputStringLen = 3*(num_atoms+num_rings+1); /*  last 3 elements are a 'zero termination' */

    if ( bCtPredecessors ) {
        if ( nOutputString = (AT_RANK *)inchi_calloc( nTotOutputStringLen, sizeof(nOutputString[0]) ) ) {
            cDelim = '-';
            for ( u = 0, k = -3 ; u < num_atoms; u ++ ) {
                k += 3;
                if ( k+6 > nTotOutputStringLen ) {
                    goto exit_error;  /* program error */
                }
                nOutputString[k]   = nNumDescendants[u]? nNumDescendants[u] : MAX_ATOMS+1;
                nOutputString[k+1] = nNum_H? 16+nNum_H[u]:0;
                nOutputString[k+2] = k? ',' : '\0';
                for ( j = 1; j <= nl[u][0] && nDfsNumber[u] > nDfsNumber[i=nl[u][j]]; j ++ ) {
                    /* closures */
                    k += 3;
                    if ( k+6 > nTotOutputStringLen ) {
                        goto exit_error;  /* program error */
                    }
                    nOutputString[k]   = i+1;  /* closure */
                    nOutputString[k+1] = 0;
                    nOutputString[k+2] = cDelim;
                }
            }
        }
    } else {
        if ( nNumDescendants ) {  /* do not need anymore */
            inchi_free( nNumDescendants );
            nNumDescendants = NULL;
        }
        /* 
            the output string contains:
              (num_atoms) atoms for the DFS (spanning) tree
              (num_atoms-1) delimiters for the DFS (spanning) tree
              1 character for each atom that has 1 terminal hydrogen atoms
              2 characters  for each atom that has 2-9 terminal hydrogen atoms
              3 characters  for each atom that has 10-99 terminal hydrogen atoms, etc.
              (num_rings) atoms for the ring closures
              (num_rings) delimiters for the ring closures
        */

        if ( nOutputString = (AT_RANK *)inchi_calloc( nTotOutputStringLen, sizeof(nOutputString[0]) ) ) {
            u               = start; /*  start atom */
            nTopStackAtom   =-1;
            memset( cNeighNumb, 0, num_atoms*sizeof(cNeighNumb[0]));
            /*  push the start atom on the stack */
            nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
            /*  output the starting atom */
            k = 0;
            nOutputString[k]   = u+1;
            nOutputString[k+1] = nNum_H? 16+nNum_H[u]:0;
            nOutputString[k+2] = '\0';

            do {
                /* advance */
                while ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i]+1,  (int)nl[i][0] >= j )
                /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
                /* replaced due to missing sequence point; undefined behavior, reported by Geoffrey Hutchison */
                {
                    k += 3;
                    if ( k+6 > nTotOutputStringLen ) {
                        goto exit_error;  /* program error */
                    }
                    cNeighNumb[i] ++;
                    u = (int)nl[i][j]; /* neighbor */

                    /* output neighbor's canonical number */
                    nOutputString[k] = u+1;

                    if ( nDfsNumber[u] > nDfsNumber[i] ) {
                        /* tree edge, 1st visit -- advance */
                        /* put 'unexplored' vertex u on the stack */
                        nStackAtom[++nTopStackAtom] = (AT_NUMB)u;

                        /* output neighbor's number of H */
                        nOutputString[k+1] = nNum_H? 16+nNum_H[u]:0;
                    } else {
                        nOutputString[k+1] = 0;
                    }
                    /* output a delimiter preceding the neighbor */
                    if ( 1 < nl[i][0] ) {
                        if ( j == 1 ) {
                            cDelim = '(';
                        } else
                        if (  j == nl[i][0] ) {
                            cDelim = ')';
                        } else {
                            cDelim = ',';
                        }
                    } else {
                        cDelim = '-';
                    }
                    nOutputString[k+2] = cDelim;
                }
                cNeighNumb[i] = 0;

                /* back up: nothing else to do */
            } while ( --nTopStackAtom >= 0 );
        }
    }
    goto exit_function;

exit_error:
    if ( nOutputString ) {
        inchi_free( nOutputString );
        nOutputString = NULL;
    }

exit_function:
    if ( nStackAtom )
        inchi_free( nStackAtom );
    if ( nNumDescendants )
        inchi_free( nNumDescendants );
    if ( nDfsNumber )
        inchi_free( nDfsNumber );
    if ( cNeighNumb )
        inchi_free( cNeighNumb );
    if ( nl )
        FreeNeighList( nl );
    return nOutputString;
}
/**********************************************************************************************/
int GetInpStructErrorType( INPUT_PARMS *ip, int err, char *pStrErrStruct, int num_inp_atoms )
{
    if ( err && err == 9 )
        return _IS_ERROR; /*  sdfile bypassed to $$$$ */
    if ( err && err < 30 )
        return _IS_FATAL;
    if ( num_inp_atoms <= 0 || err ) {
        if ( 98 == err && 0 == num_inp_atoms && ip->bAllowEmptyStructure )
            return _IS_WARNING;
        return _IS_ERROR;
    }
    if ( pStrErrStruct[0] )
        return _IS_WARNING;
    return _IS_OKAY;
}
/**********************************************************************************************/
int ProcessStructError( INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *log_file, /*int err,*/ 
                        char *pStrErrStruct, int nErrorType,
                        int *bXmlStructStarted, long num_inp, INPUT_PARMS *ip, char *pStr, int nStrLen )
{
    int b_ok;
#ifdef INCHI_LIB
    int bPlainText = (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) &&
                     (ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
                    !(ip->bINChIOutputOptions & INCHI_OUT_XML);
#else
    int bPlainText = 0;
#endif
    if ( !bPlainText && *bXmlStructStarted <= 0 ) {
        return nErrorType;
    }
    /*  Fatal error, Error, Warning */
    if ( nErrorType ) {
        if ( bPlainText ) {
            if ( !(b_ok=OutputINChIPlainError( output_file, pStr, nStrLen, pStrErrStruct, nErrorType ) ) ) {
                inchi_ios_eprint( log_file, "Cannot create message for error (structure #%ld.%s%s%s%s) Terminating.\n",
                                                            num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            } else {
                inchi_ios_print( output_file, "\n" ); /* add a blank line after the WINCHI Window message */
            }
        } else {
            if ( !(b_ok=OutputINChIXmlError( output_file, pStr, nStrLen, 2, /*err,*/ pStrErrStruct, nErrorType ) ) ) {
                inchi_ios_eprint( log_file, "Cannot create xml tag for error (structure #%ld.%s%s%s%s) Terminating.\n",
                                                            num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            }
            if ( !b_ok || nErrorType == _IS_FATAL || nErrorType == _IS_ERROR ) {
                /*  close current structure output */
                if ( !OutputINChIXmlStructEndTag( output_file, pStr, nStrLen, 1 ) ) {
                    inchi_ios_eprint( log_file, "Cannot create end xml tag for structure #%ld.%s%s%s%s Terminating.\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                    *bXmlStructStarted = -1;
                    b_ok = 0;
                } else {
                    *bXmlStructStarted = 0;
                }
            }
        }
        return b_ok? nErrorType : _IS_FATAL;
    }
    return nErrorType;
    
}

#if( TEST_RENUMB_ATOMS == 1 ) /*  { */
/***************************************************************************************/
int CompareStereoINChI( INChI_Stereo *s1, INChI_Stereo *s2 )
{
    if ( s1 == NULL && s2 == NULL )
        return 0;
    if ( (s1 == NULL) ^ (s2 == NULL) )
        return 20;

    if ( s1->nNumberOfStereoCenters != s2->nNumberOfStereoCenters )
        return 21;
    if ( s1->nNumberOfStereoCenters > 0 ) {
        if ( memcmp( s1->nNumber, s2->nNumber, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
            return 22;
        if ( memcmp( s1->t_parity, s2->t_parity, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
            return 23;
        if ( s1->nNumberInv && s2->nNumberInv ) {
            if ( memcmp( s1->nNumberInv, s2->nNumberInv, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
                return 28;
            if ( memcmp( s1->t_parityInv, s2->t_parityInv, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
                return 29;
            if ( s1->nCompInv2Abs != s2->nCompInv2Abs ||
                 s1->bTrivialInv  != s2->bTrivialInv ) {
                return 30;
            }
        } else
        if ( s1->nNumberInv || s2->nNumberInv ) {
            return 31;
        }
    }
    if ( s1->nNumberOfStereoBonds != s2->nNumberOfStereoBonds )
        return 24;
    if ( s1->nNumberOfStereoBonds > 0 ) {
        if ( memcmp( s1->nBondAtom1, s2->nBondAtom1, s1->nNumberOfStereoBonds*sizeof(s1->nBondAtom1[0]) ) )
            return 25;
        if ( memcmp( s1->nBondAtom2, s2->nBondAtom2, s1->nNumberOfStereoBonds*sizeof(s1->nBondAtom2[0]) ) )
            return 26;
        if ( memcmp( s1->b_parity, s2->b_parity, s1->nNumberOfStereoBonds*sizeof(s1->b_parity[0]) ) )
            return 27;
    }
    return 0;
}
/***************************************************************************************/
int CompareINChI( INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2 )
{
    int ret;
    if ( i1 == NULL && i2 == NULL )
        return 0;
    if ( (i1 == NULL) ^ (i2 == NULL) )
        return 1;
    
    if ( i1->nErrorCode == i2->nErrorCode ) {
        if ( i1->nErrorCode )
            return 0;
    } else {
        return 2;
    }
    
    if ( i1->nNumberOfAtoms != i2->nNumberOfAtoms )
        return 3;
    if ( i1->nNumberOfAtoms > 0 ) {
        if ( memcmp( i1->nAtom, i2->nAtom, i1->nNumberOfAtoms*sizeof(i1->nAtom[0]) ) )
            return 4;
        if ( memcmp( i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms*sizeof(i1->nNum_H[0]) ) )
            return 5;
        if ( i1->nNum_H_fixed && i2->nNum_H_fixed && 
             memcmp( i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms*sizeof(i1->nNum_H_fixed[0]) ) ) {
            return 6;
        }
        if ( strcmp( i1->szHillFormula, i2->szHillFormula ) )
            return 7;
    }

    if ( i1->lenConnTable != i2->lenConnTable )
        return 8;
    if ( i1->lenConnTable > 0 && memcmp( i1->nConnTable, i2->nConnTable, i1->lenConnTable*sizeof(i1->nConnTable[0]) ) )
        return 9;

    if ( i1->lenTautomer != i2->lenTautomer )
        return 10;
    if ( i1->lenTautomer > 0 && memcmp( i1->nTautomer, i2->nTautomer, i1->lenTautomer*sizeof(i1->nTautomer[0]) ) )
        return 11;

    if ( i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms )
        return 12;
    if ( i1->nNumberOfIsotopicAtoms > 0 && memcmp( i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms*sizeof(i1->IsotopicAtom[0]) ) )
        return 13;

    if ( i1->nNumberOfIsotopicTGroups != i2->nNumberOfIsotopicTGroups )
        return 14;
    if ( i1->nNumberOfIsotopicTGroups > 0 && memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup, i1->nNumberOfIsotopicTGroups*sizeof(i1->IsotopicTGroup[0]) ) )
        return 15;
    if ( a1->nNumRemovedProtons != a2->nNumRemovedProtons )
        return 16;
    if ( memcmp( a1->nNumRemovedIsotopicH, a2->nNumRemovedIsotopicH, sizeof(a1->nNumRemovedIsotopicH) ) )
        return 17;
    if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
        if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
             memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                    sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
            return 18;
    } else
    if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
        return 19;
    }
    /* ret = 20..31 */
    if ( ret = CompareStereoINChI( i1->Stereo, i2->Stereo ) )
        return ret;
    /* ret = 40..51 */
    if ( ret = CompareStereoINChI( i1->StereoIsotopic, i2->StereoIsotopic ) )
        return ret+20;

    return 0;
}
#endif /*  } TEST_RENUMB_ATOMS == 1 */
#if( READ_INCHI_STRING == 1 ) /*  { */
/*************************************************************************************/
int CompareReversedStereoINChI( INChI_Stereo *s1/* InChI from reversed struct */, INChI_Stereo *s2 /* input InChI */)
{
    if ( s1 == NULL && s2 == NULL )
        return 0;
    if ( (s1 == NULL) ^ (s2 == NULL) ) {
        INChI_Stereo *s = s1? s1 : s2;
        if ( s->nNumberOfStereoCenters || s->nNumberOfStereoBonds ) {
            return 20; /* Diff: Missing Stereo */
        } else {
            return 0;
        }
    }

    if ( s1->nNumberOfStereoCenters != s2->nNumberOfStereoCenters )
        return 21;      /* Diff: Number of sp3 stereocenters */
    if ( s1->nNumberOfStereoCenters > 0 ) {
        if ( memcmp( s1->nNumber, s2->nNumber, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
            return 22;  /* Diff: sp3 stereocenter locations */
        if ( memcmp( s1->t_parity, s2->t_parity, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
            return 23;  /* Diff: sp3 stereocenter parities */
        if ( s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs )
            return 24;  /* Diff: sp3 inversion */
        /*
        if ( s1->nNumberInv && s2->nNumberInv ) {
            if ( memcmp( s1->nNumberInv, s2->nNumberInv, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
                return 25;
            if ( memcmp( s1->t_parityInv, s2->t_parityInv, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
                return 26;
            if ( s1->nCompInv2Abs != s2->nCompInv2Abs ||
                 s1->bTrivialInv  != s2->bTrivialInv ) {
                return 27;
            }
        } else
        if ( s1->nNumberInv || s2->nNumberInv ) {
            return 28;
        }
        */
    }
    if ( s1->nNumberOfStereoBonds != s2->nNumberOfStereoBonds )
        return 25;      /* Diff: Number of stereobonds */
    if ( s1->nNumberOfStereoBonds > 0 ) {
        if ( memcmp( s1->nBondAtom1, s2->nBondAtom1, s1->nNumberOfStereoBonds*sizeof(s1->nBondAtom1[0]) ) )
            return 26; /* Diff: Stereobond 1st atom locations */
        if ( memcmp( s1->nBondAtom2, s2->nBondAtom2, s1->nNumberOfStereoBonds*sizeof(s1->nBondAtom2[0]) ) )
            return 27; /* Diff: Stereobond 2nd atom locations */
        if ( memcmp( s1->b_parity, s2->b_parity, s1->nNumberOfStereoBonds*sizeof(s1->b_parity[0]) ) )
            return 28; /* Diff: Stereobond parities */
    }
    return 0;
}
/*************************************************************************************/
int CompareReversedStereoINChI2( INChI_Stereo *s1/* InChI from reversed struct */, INChI_Stereo *s2 /* input InChI */, ICR *picr)
{
    int ret = 0;
    int j1, j2, num_eq, num_dif, num_extra_undf, num_miss_undf, num_in1_only, num_in2_only;
    int bAddSb = !(picr->num_sb_undef_in1_only + picr->num_sb_in1_only + picr->num_sb_in2_only);
    int bAddSc = !(picr->num_sc_undef_in1_only + picr->num_sc_in1_only + picr->num_sc_in2_only);
    
    int nNumSc1 = s1? s1->nNumberOfStereoCenters : 0;
    int nNumSc2 = s2? s2->nNumberOfStereoCenters : 0;
    int nNumSb1 = s1? s1->nNumberOfStereoBonds   : 0;
    int nNumSb2 = s2? s2->nNumberOfStereoBonds   : 0;
    
    if ( (nNumSc1 || nNumSc1) &&
         ( nNumSc1 != nNumSc2 ||
           memcmp( s1->nNumber,  s2->nNumber,  nNumSc1*sizeof(s1->nNumber[0] ) ) ||
           memcmp( s1->t_parity, s2->t_parity, nNumSc1*sizeof(s1->t_parity[0]) ) ) ) {

        num_eq = num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < nNumSc1 && j2 < nNumSc2; ) {
            if ( s1->nNumber[j1] ==  s2->nNumber[j2] ) {
                if ( s1->t_parity[j1] == s2->t_parity[j2] ) {
                    num_eq ++;
                } else {
                    num_dif ++;
                }
                j1 ++;
                j2 ++;
            } else
            if ( s1->nNumber[j1] < s2->nNumber[j2] ) {
                num_in1_only ++;
                if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                    num_extra_undf ++;
                }
                if ( bAddSc ) {
                    if ( picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY )
                        picr->sc_in1_only[picr->num_sc_in1_only ++] = j1;
                    if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                        if ( picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF )
                            picr->sc_undef_in1_only[picr->num_sc_undef_in1_only ++] = j1;
                    }
                }
                j1 ++;
            } else {
                num_in2_only ++;
                if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                    num_miss_undf ++;
                }
                if ( bAddSc ) {
                    if ( picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY )
                        picr->sc_in2_only[picr->num_sc_in2_only ++] = j2;
                    if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                        if ( picr->num_sc_undef_in2_only < ICR_MAX_SC_UNDF )
                            picr->sc_undef_in2_only[picr->num_sc_undef_in2_only ++] = j1;
                    }
                }
                j2 ++;
            }
        }
        while ( j1 < nNumSc1 ) {
            if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                num_extra_undf ++;
            }
            num_in1_only ++;
            if ( bAddSc ) {
                if ( picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY )
                    picr->sc_in1_only[picr->num_sc_in1_only ++] = j1;
                if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                    if ( picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF )
                        picr->sc_undef_in1_only[picr->num_sc_undef_in1_only ++] = j1;
                }
            }
            j1 ++;
        }
        while ( j2 < nNumSc2 ) {
            if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                num_miss_undf ++;
            }
            num_in2_only ++;
            if ( bAddSc ) {
                if ( picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY )
                    picr->sc_in2_only[picr->num_sc_in2_only ++] = j2;
            }
            j2 ++;
        }
        if ( num_dif ) {
            ret |= IDIF_SC_PARITY; 
        }
        if ( num_in1_only ) {
            if ( num_extra_undf ) {
                ret |= IDIF_SC_EXTRA_UNDF;
            }
            if ( num_in1_only != num_extra_undf ) {
                ret |= IDIF_SC_EXTRA;
            }
        }
        if ( num_in2_only ) {
            if ( num_miss_undf ) {
                ret |= IDIF_SC_MISS_UNDF;
            }
            if ( num_in2_only != num_miss_undf ) {
                ret |= IDIF_SC_MISS;
            }
        }
    }
    if ( s1 && s2 && s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs ) {
        ret |= IDIF_SC_INV;
    }

    if ( (nNumSb1 || nNumSb2 ) &&
         (nNumSb1 != nNumSb2 ||
          memcmp( s1->nBondAtom1, s2->nBondAtom1, nNumSb1*sizeof(s1->nBondAtom1[0]) ) ||
          memcmp( s1->nBondAtom2, s2->nBondAtom2, nNumSb1*sizeof(s1->nBondAtom2[0]) ) ||
          memcmp( s1->b_parity,   s2->b_parity,   nNumSb1*sizeof(s1->b_parity[0]) ) ) ) {

        num_eq = num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < nNumSb1 && j2 < nNumSb2; ) {
            if ( s1->nBondAtom1[j1] ==  s2->nBondAtom1[j2] &&
                 s1->nBondAtom2[j1] ==  s2->nBondAtom2[j2] ) {
                if ( s1->b_parity[j1] == s2->b_parity[j2] ) {
                    num_eq ++;
                } else {
                    num_dif ++;
                }
                j1 ++;
                j2 ++;
            } else
            if ( s1->nBondAtom1[j1] <  s2->nBondAtom1[j2] ||
                 s1->nBondAtom1[j1] == s2->nBondAtom1[j2] && s1->nBondAtom2[j1] <  s2->nBondAtom2[j2]) {
                num_in1_only ++;
                if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                    num_extra_undf ++;
                }
                if ( bAddSb ) {
                    if ( picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY )
                        picr->sb_in1_only[picr->num_sb_in1_only ++] = j1;
                    if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                        if ( picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF )
                            picr->sb_undef_in1_only[picr->num_sb_undef_in1_only ++] = j1;
                    }
                }
                j1 ++;
            } else {
                num_in2_only ++;
                if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                    num_miss_undf ++;
                }
                if ( bAddSb ) {
                    if ( picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY )
                        picr->sb_in2_only[picr->num_sb_in2_only ++] = j2;
                    if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                        if ( picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF )
                            picr->sb_undef_in2_only[picr->num_sb_undef_in2_only ++] = j1;
                    }
                }
                j2 ++;
            }
        }
        while ( j1 < nNumSb1 ) {
            num_in1_only ++;
            if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                num_extra_undf ++;
            }
            if ( bAddSb ) {
                if ( picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY )
                    picr->sb_in1_only[picr->num_sb_in1_only ++] = j1;
                if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                    if ( picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF )
                        picr->sb_undef_in1_only[picr->num_sb_undef_in1_only ++] = j1;
                }
            }
            j1 ++;
        }
        while ( j2 < nNumSb2 ) {
            num_in2_only ++;
            if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                num_miss_undf ++;
            }
            if ( bAddSb ) {
                if ( picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY )
                    picr->sb_in2_only[picr->num_sb_in2_only ++] = j2;
                if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                    if ( picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF )
                        picr->sb_undef_in2_only[picr->num_sb_undef_in2_only ++] = j1;
                }
            }
            j2 ++;
        }
        if ( num_dif ) {
            ret |= IDIF_SB_PARITY; 
        }
        if ( num_in1_only ) {
            if ( num_extra_undf ) {
                ret |= IDIF_SB_EXTRA_UNDF;
            }
            if ( num_in1_only != num_extra_undf ) {
                ret |= IDIF_SB_EXTRA;
            }
        }
        if ( num_in2_only ) {
            if ( num_miss_undf ) {
                ret |= IDIF_SB_MISS_UNDF;
            }
            if ( num_in2_only != num_miss_undf ) {
                ret |= IDIF_SB_MISS;
            }
        }
    }

    return ret;
}
/*************************************************************************************/
int CompareReversedINChI( INChI *i1 /* InChI from reversed struct */, INChI *i2 /* input InChI */, INChI_Aux *a1, INChI_Aux *a2 )
{
    int ret;
    if ( i1 == NULL && i2 == NULL )
        return 0;
    if ( (i1 == NULL) ^ (i2 == NULL) )
        return 1; /* Diff: Missing InChI */
    
    if ( i1->nErrorCode == i2->nErrorCode ) {
        if ( i1->nErrorCode )
            return 0;
    } else {
        return 2; /* Diff: Error codes */
    }
    if ( i1->bDeleted != i2->bDeleted ) {
        return 1; /* Diff: Missing InChI */
    }
    if ( i1->nNumberOfAtoms != i2->nNumberOfAtoms )
        return 3;  /* Diff: Num. atoms */
    if ( i1->nNumberOfAtoms > 0 ) {
        if ( memcmp( i1->nAtom, i2->nAtom, i1->nNumberOfAtoms*sizeof(i1->nAtom[0]) ) )
            return 4; /* Diff: Elements */
        if ( strcmp( i1->szHillFormula, i2->szHillFormula ) )
            return 7; /* Diff: Hill Formulas */
        if ( memcmp( i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms*sizeof(i1->nNum_H[0]) ) ) {
            if ( i1->lenConnTable > 1 || i2->lenConnTable > 1 ) {
                return 5; /* Diff: H Locations (mobile H present) */
            } else {
                return 6; /* Diff: H Locations (no mobile H) */
            }
        }
        /* fixed H */
        if ( i1->nNum_H_fixed || i2->nNum_H_fixed ) {
            int bHasFixedH1 = 0, bHasFixedH2 = 0, i, j1, j2;
            if ( i1->nNum_H_fixed ) {
                for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] ) {
                        bHasFixedH1 ++;
                    }
                }
            }
            if ( i2->nNum_H_fixed ) {
                for ( i = 0; i < i2->nNumberOfAtoms; i ++ ) {
                    if ( i2->nNum_H_fixed[i] ) {
                        bHasFixedH2 ++;
                    }
                }
            }
            /* count the differences */
            j1 = j2 = 0;
            if ( bHasFixedH1 && !bHasFixedH2 ) {
                for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] > 0 ) {
                        j1 ++;
                    } else
                    if ( i1->nNum_H_fixed[i] < 0 ) {
                        j2 ++;
                    }
                }

                return 18; /* Diff: Extra Fixed-H */
            } else
            if ( !bHasFixedH1 && bHasFixedH2 ) {
                for ( i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( 0 > i2->nNum_H_fixed[i] ) {
                        j1 ++;
                    } else
                    if ( 0 < i2->nNum_H_fixed[i] ) {
                        j2 ++;
                    }
                }
                return 19; /* Diff: Missed Fixed-H */
            } else
            if ( bHasFixedH1 && bHasFixedH2 &&
                 memcmp( i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms*sizeof(i1->nNum_H_fixed[0]) ) ) {
                for ( i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i] ) {
                        j1 ++;
                    } else
                    if ( i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i] ) {
                        j2 ++;
                    }
                }
            }
            ret = (j1 && j2)? 20 : j1? 18 : j2? 19 : 0;
            if ( ret ) {
                return ret; /* 20 => Diff: NotEql Fixed-H */
                            /* 19 => Diff: Missed Fixed-H (i1 has less) */
                            /* 18 => Diff: Extra Fixed-H  (i1 has more) */
            }
        }
    }

    if ( i1->lenConnTable != i2->lenConnTable )
        return 8; /* Diff: Connections length */
    if ( i1->lenConnTable > 0 && memcmp( i1->nConnTable, i2->nConnTable, i1->lenConnTable*sizeof(i1->nConnTable[0]) ) )
        return 9; /* Diff: Connections */
    /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    if ( i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1) )
        return 10; /* Diff: Mobile groups length */ /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
    if ( (i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
         memcmp( i1->nTautomer, i2->nTautomer, i1->lenTautomer*sizeof(i1->nTautomer[0]) ) )
        return 11; /* Diff: Mobile groups */

    if ( i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms )
        return 12; /* Diff: Isotopic atoms number */
    if ( i1->nNumberOfIsotopicAtoms > 0 && memcmp( i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms*sizeof(i1->IsotopicAtom[0]) ) )
        return 13; /* Diff: Isotopic atoms */
    if ( i1->nTotalCharge != i2->nTotalCharge )
        return 14; /* Diff: Charge */
/*
    if ( i1->nNumberOfIsotopicTGroups != i2->nNumberOfIsotopicTGroups )
        return 14;
    if ( i1->nNumberOfIsotopicTGroups > 0 && memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup, i1->nNumberOfIsotopicTGroups*sizeof(i1->IsotopicTGroup[0]) ) )
        return 15;
*/
    if ( a1 && a2 ) {
        if ( a1->nNumRemovedProtons != a2->nNumRemovedProtons )
            return 16; /* Diff: Number of removed protons */
        if ( memcmp( a1->nNumRemovedIsotopicH, a2->nNumRemovedIsotopicH, sizeof(a1->nNumRemovedIsotopicH) ) )
            return 17; /* Diff: Removed isotopic H */
    }
/*
    if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
        if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
             memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                    sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
            return 18;
    } else
    if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
        return 19;
    }
*/
    /* ret = 20..31 => 40..51 */
    if ( ret = CompareReversedStereoINChI( i1->Stereo, i2->Stereo ) )
        return ret+20;
    /* ret = 40..51 => 60..71 */

    if ( !i2->StereoIsotopic && i2->Stereo && i1->StereoIsotopic &&
         0 < (i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters) &&
         0 == CompareReversedStereoINChI( i1->StereoIsotopic, i2->Stereo ) ) {
        /* InChI from reversed structure does not contain fully duplicated isotopic stereo */
        ;
    } else

    if ( ret = CompareReversedStereoINChI( i1->StereoIsotopic, i2->StereoIsotopic ) ) {
        return ret+40;
    }

    return 0;
}

/*******************************************************************************/
int CompareIcr( ICR *picr1, ICR *picr2, INCHI_MODE *pin1, INCHI_MODE *pin2, INCHI_MODE mask )
{
    int nNumExtraBits1 = 0, nNumExtraBits2 = 0, bit1, bit2;
    INCHI_MODE Flg1=picr1->flags, Flg2 = picr2->flags, cur_bit = 1, in1, in2;
    int i, ret;

    /* compare flags */
    in1 = in2 = 0;
    for ( i = 0; Flg1 || Flg2; i ++, Flg1 >>= 1, Flg2 >>= 1, cur_bit <<= 1 ) {
        if ( !(mask & cur_bit) ) {
            continue;
        }
        bit1 = Flg1 & 1;
        bit2 = Flg2 & 1;
        if ( bit1 && !bit2 ) {
            in1 |= 1 << i;
            nNumExtraBits1 ++;
        } else
        if ( !bit1 && bit2 ) {
            in2 |= 1 << i;
            nNumExtraBits2 ++;
        }
    }
    if ( nNumExtraBits1 && !nNumExtraBits2 ) {
        ret = 1;
    } else
    if ( !nNumExtraBits1 && nNumExtraBits2 ) {
        ret = -1;
    } else
    if ( !in1 && !in2 ) {
        ret = 0;
    } else {
        ret = 2; /* compare produced undefined results */
    }
    if ( pin1 ) *pin1 = in1;
    if ( pin2 ) *pin2 = in2;
    /* more detailed compare not implemented */
    return ret;
}

/*********************************************************************************************************/
INCHI_MODE CompareReversedINChI2( INChI *i1 /* InChI from reversed struct */, INChI *i2 /* input InChI */,
                                  INChI_Aux *a1, INChI_Aux *a2, ICR *picr, int *err )
{
    INCHI_MODE ret = 0;
    INChI_Stereo *Stereo1=NULL, *Stereo2=NULL;
    int  n1, n2, m, j, j1, j2, ret2, num_H1, num_H2;
    
    *err = 0;

    memset( picr, 0, sizeof(*picr) );

    if ( i1 == NULL && i2 == NULL )
        return 0;
    if ( (i1 == NULL) ^ (i2 == NULL) ) {
        ret |= IDIF_PROBLEM; /* one InChI exists while another doesn't */
        goto exit_function;
    }
    
    if ( i1->nErrorCode == i2->nErrorCode ) {
        if ( i1->nErrorCode ) {
            ret |= IDIF_PROBLEM; /* both InChI have same error codes */
            goto exit_function;
        }
    } else {
        ret |= IDIF_PROBLEM; /* at least one InChI has an error code */
        goto exit_function;
    }
    
    if ( i1->nNumberOfAtoms != i2->nNumberOfAtoms ) {
        ret |= IDIF_NUM_AT;
        goto exit_function;
    }
    if ( i1->nNumberOfAtoms > 0 ) {
        if ( memcmp( i1->nAtom, i2->nAtom, i1->nNumberOfAtoms*sizeof(i1->nAtom[0]) ) ) {
            ret |= IDIF_ATOMS;
            goto exit_function;
        }
        /* IDIF_NON_TAUT_H,  IDIF_MORE_FH, IDIF_LESS_FH */
        if ( memcmp( i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms*sizeof(i1->nNum_H[0]) ) ) {
            ret |= IDIF_POSITION_H;
            for ( j1 = 0; j1 < i1->nNumberOfAtoms; j1 ++ ) {
                if ( i1->nNum_H[j1] != i2->nNum_H[j1] && picr->num_diff_pos_H < ICR_MAX_DIFF_FIXED_H ) {
                    picr->diff_pos_H_at[picr->num_diff_pos_H] = j1;
                    picr->diff_pos_H_nH[picr->num_diff_pos_H] = i1->nNum_H[j1] - i2->nNum_H[j1];
                    picr->num_diff_pos_H ++;
                }
            }
        }
        /* fixed H */
        if ( i1->nNum_H_fixed || i2->nNum_H_fixed ) {
            int bHasFixedH1 = 0, bHasFixedH2 = 0, i;
            if ( i1->nNum_H_fixed ) {
                for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] ) {
                        bHasFixedH1 ++;
                    }
                }
            }
            if ( i2->nNum_H_fixed ) {
                for ( i = 0; i < i2->nNumberOfAtoms; i ++ ) {
                    if ( i2->nNum_H_fixed[i] ) {
                        bHasFixedH2 ++;
                    }
                }
            }
            if ( bHasFixedH1 && !bHasFixedH2 ) {
                for ( i = j = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] ) {
                        if ( j < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at1_more[j] = i;
                            picr->fixed_H_nH1_more[j] = i1->nNum_H_fixed[i];
                            j ++;
                        }
                    }
                }
                picr->num_fixed_H1_more = j;
                ret |= IDIF_MORE_FH; /* Extra Fixed-H */
            } else
            if ( !bHasFixedH1 && bHasFixedH2 ) {
                for ( i = j = 0; i < i2->nNumberOfAtoms; i ++ ) {
                    if ( i2->nNum_H_fixed[i] ) {
                        if ( j < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at2_more[j] = i;
                            picr->fixed_H_nH2_more[j] = i2->nNum_H_fixed[i];
                            j ++;
                        }
                    }
                }
                picr->num_fixed_H2_more = j;
                ret |= IDIF_LESS_FH; /* Missed Fixed-H */
            } else
            if ( bHasFixedH1 && bHasFixedH2 &&
                 memcmp( i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms*sizeof(i1->nNum_H_fixed[0]) ) ) {
                for ( i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i] ) {
                        if ( j1 < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at1_more[j1] = i;
                            picr->fixed_H_nH1_more[j1] = i1->nNum_H_fixed[i] - i2->nNum_H_fixed[i];
                            j1 ++;
                        }
                    } else
                    if ( i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i] ) {
                        if ( j2 < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at2_more[j2] = i;
                            picr->fixed_H_nH2_more[j2] = i2->nNum_H_fixed[i] - i1->nNum_H_fixed[i];
                            j2 ++;
                        }
                    }
                }
                ret |= (j1? IDIF_MORE_FH:0) | (j2? IDIF_LESS_FH:0);
                picr->num_fixed_H1_more = j1;
                picr->num_fixed_H2_more = j2;
            }
        }
    }
    /* compare formulas and H */
    num_H1 = 0;
    num_H2 = 0;
    ret2 = CompareHillFormulasNoH( i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2 );
    picr->tot_num_H1 = num_H1;
    picr->tot_num_H2 = num_H2;
    if ( ret2 ) {
        ret |= IDIF_NUM_EL;
        goto exit_function;
    }
    if ( num_H1 > num_H2 ) {
        ret |= IDIF_MORE_H;
    }
    if ( num_H1 < num_H2 ) {
        ret |= IDIF_LESS_H;
    }

    if ( i1->lenConnTable != i2->lenConnTable ) {
        ret |= IDIF_CON_LEN;
        goto exit_function;
    } else
    if ( i1->lenConnTable > 0 && memcmp( i1->nConnTable, i2->nConnTable, i1->lenConnTable*sizeof(i1->nConnTable[0]) ) ) {
        ret |= IDIF_CON_TBL;
        goto exit_function;
    }
    /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
/*
    if ( i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1) ) {
        ret |=  IDIF_TAUT_LEN; 
    }
*/
    /* compare number of t-groups */
    n1 = i1->lenTautomer? i1->nTautomer[0] : 0;
    n2 = i2->lenTautomer? i2->nTautomer[0] : 0;
    if ( !n1 && n2 ) {
        ret |= IDIF_NO_TAUT;
    } else
    if ( n1 && !n2 ) {
        ret |= IDIF_WRONG_TAUT;
    } else
    if ( n1 == 1 && n2 > 1 ) {
        ret |= IDIF_SINGLE_TG;
    } else
    if ( n1 > 1 && n2 == 1 ) {
        ret |= IDIF_MULTIPLE_TG;
    } else
    if ( n1 != n2 ) {
        ret |= IDIF_NUM_TG;
    }
    if ( n1 || n2 ) {
        /* number of endpoints */
        int num1 = 0, num2 = 0, num_M1=0, num_M2=0;
        int len, num_eq, num_in1_only, num_in2_only;
        AT_NUMB *pe1 = (AT_NUMB *)inchi_malloc( (i1->lenTautomer+1) * sizeof(pe1[0]) );
        AT_NUMB *pe2 = (AT_NUMB *)inchi_malloc( (i2->lenTautomer+1) * sizeof(pe2[0]) );
        num_H1 = num_H2=0;
        /* collect endpoints, H, (-) */
        if ( !pe1 || !pe2 ) {
            if ( pe1 ) inchi_free( pe1 );
            if ( pe2 ) inchi_free( pe2 );
            *err = -1; /* allocation error */
            goto exit_function;
        }
        for ( m = 1; m < i1->lenTautomer; m += len ) {
            len = i1->nTautomer[m ++];
            num_H1 += i1->nTautomer[m];
            num_M1 += i1->nTautomer[m+1];
            for ( j = 2; j < len; j ++ ) {
                pe1[num1 ++] = i1->nTautomer[m + j];
            }
        }
        for ( m = 1; m < i2->lenTautomer; m += len ) {
            len = i2->nTautomer[m ++];
            num_H2 += i2->nTautomer[m];
            num_M2 += i2->nTautomer[m+1];
            for ( j = 2; j < len; j ++ ) {
                pe2[num2 ++] = i2->nTautomer[m + j];
            }
        }
        picr->num_taut_H1 = num_H1;
        picr->num_taut_H2 = num_H2;
        picr->num_taut_M1 = num_M1;
        picr->num_taut_M2 = num_M2;
        /* sort endpoints */
        insertions_sort_AT_RANK( pe1, num1 );
        insertions_sort_AT_RANK( pe2, num2 );
        /* compare */
        /*
        if ( num1 < num2 ) {
            ret |= IDIF_LESS_TG_ENDP;
        } else
        if ( num1 > num2 ) {
            ret |= IDIF_MORE_TG_ENDP;
        }
        */
        /* compare all */
        num_eq = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < num1 && j2 < num2; ) {
            if( pe1[j1] == pe2[j2] ) {
                j1 ++;
                j2 ++;
                num_eq ++;
            } else
            if ( pe1[j1] < pe2[j2] ) { /* BC: fixed, was pe2[j1] 2006-03-27 */
                if ( picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY ) {
                    picr->endp_in1_only[picr->num_endp_in1_only ++] = pe1[j1];
                }
                j1 ++;
                num_in1_only ++;
            } else {
                if ( picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY ) {
                    picr->endp_in2_only[picr->num_endp_in2_only ++] = pe2[j2];
                }
                j2 ++;
                num_in2_only ++;
            }
        }
        while ( j1 < num1 ) {
            if ( picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY ) {
                picr->endp_in1_only[picr->num_endp_in1_only ++] = pe1[j1];
            }
            j1 ++;
            num_in1_only ++;
        }
        while ( j2 < num2 ) {
            if ( picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY ) {
                picr->endp_in2_only[picr->num_endp_in2_only ++] = pe2[j2];
            }
            j2 ++;
            num_in2_only ++;
        }
        if ( num_in1_only ) {
            ret |= IDIF_EXTRA_TG_ENDP;
        }
        if ( num_in2_only ) {
            ret |= IDIF_MISS_TG_ENDP;
        }
        if ( !num_in1_only && !num_in2_only && num_eq ) {
           ; /* same t-groups endpoints */
        } else {
           ret |= IDIF_DIFF_TG_ENDP;
        }
        inchi_free( pe1 );
        inchi_free( pe2 );

    }

    if ( (i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
         ( i1->lenTautomer != i2->lenTautomer ||
         memcmp( i1->nTautomer, i2->nTautomer, i1->lenTautomer*sizeof(i1->nTautomer[0]) ) ) )
        ret |= IDIF_TG;

    if ( i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms ) {
        ret |= IDIF_NUM_ISO_AT;
    } else
    if ( i1->nNumberOfIsotopicAtoms > 0 && memcmp( i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms*sizeof(i1->IsotopicAtom[0]) ) )
        ret |= IDIF_ISO_AT;
    if ( i1->nTotalCharge != i2->nTotalCharge )
        ret |= IDIF_CHARGE;
    if ( a1 && a1->nNumRemovedProtons && (!a2 || a2->nNumRemovedProtons != a1->nNumRemovedProtons) ) {
        ret |= IDIF_REM_PROT;
    }
    if ( a1 && (!a2 || 
         a2->nNumRemovedIsotopicH[0] != a1->nNumRemovedIsotopicH[0] ||
         a2->nNumRemovedIsotopicH[1] != a1->nNumRemovedIsotopicH[1] ||
         a2->nNumRemovedIsotopicH[2] != a1->nNumRemovedIsotopicH[2]) ) {
        ret |= IDIF_REM_ISO_H;
    }

/*
    if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
        if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
             memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                    sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
            return 18;
    } else
    if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
        return 19;
    }
*/
    if ( i1->StereoIsotopic &&
         i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters ) {
        Stereo1 = i1->StereoIsotopic;
    } else {
        Stereo1 = i1->Stereo;
    }
    if ( i2->StereoIsotopic &&
         i2->StereoIsotopic->nNumberOfStereoBonds + i2->StereoIsotopic->nNumberOfStereoCenters ) {
        Stereo2 = i2->StereoIsotopic;
    } else {
        Stereo2 = i2->Stereo;
    }
    ret |= CompareReversedStereoINChI2( Stereo1, Stereo2, picr );

exit_function:

    picr->flags = ret;

    return ret;
}
#endif  /* } READ_INCHI_STRING */
/***************************************************************************************/
int  Create_INChI( INChI **ppINChI, INChI_Aux **ppINChI_Aux, ORIG_ATOM_DATA *orig_inp_data, /* not used */
                  inp_ATOM *inp_at, INP_ATOM_DATA *out_norm_data[2],
                  int num_inp_at, INCHI_MODE nUserMode,
                  INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone,
                  struct tagInchiTime *ulMaxTime, T_GROUP_INFO *ti_out, char *pStrErrStruct)
{
/*
#define NON_TAUT 0
#define TAUT     1    
*/
    sp_ATOM  *at[TAUT_NUM]; /* at[0]=>non-tautomeric, at[1]=>tautomeric */
    /* inp_ATOM *out_norm_taut_at, *out_norm_nontaut_at; */
    int                       i, n1, n2, num_atoms, num_at_tg, num_removed_H, num_removed_H_taut=0, ret=0, ret2=0;
    INCHI_MODE                 nMode=0;
    T_GROUP_INFO              vt_group_info;
    T_GROUP_INFO              vt_group_info_orig;
    T_GROUP_INFO * /*const*/  t_group_info        = &vt_group_info;                
    T_GROUP_INFO * /*const*/  t_group_info_orig   = &vt_group_info_orig;
    
    CANON_STAT  CS, CS2;
    CANON_STAT *pCS  = &CS;
    CANON_STAT *pCS2 = &CS2;  /*  save all allocations to avoid memory leaks in case Canon_INChI() removes the pointer */
    
    ATOM_SIZES  s[TAUT_NUM];

    BCN Bcn;
    BCN *pBCN = &Bcn;

    int bHasIsotopicAtoms  = 0;
    int bMayHaveStereo     = 0;
    int num_taut_at        = 0;

    inp_ATOM *out_at = NULL; /*, *norm_at_fixed_bonds[TAUT_NUM]; */ /*  = {out_norm_nontaut_at, out_norm_taut_at} ; */
    INChI     *pINChI=NULL;      /* added initialization 2006-03 */
    INChI_Aux *pINChI_Aux=NULL;  /* added initialization 2006-03 */
    int        bPointedEdgeStereo = ((TG_FLAG_POINTED_EDGE_STEREO & *pbTautFlags)? PES_BIT_POINT_EDGE_STEREO:0)
                                  | ((TG_FLAG_PHOSPHINE_STEREO    & *pbTautFlags)? PES_BIT_PHOSPHINE_STEREO :0)
                                  | ((TG_FLAG_ARSINE_STEREO       & *pbTautFlags)? PES_BIT_ARSINE_STEREO    :0)
                                  | ((TG_FLAG_FIX_SP3_BUG         & *pbTautFlags)? PES_BIT_FIX_SP3_BUG      :0);
    INCHI_MODE bTautFlags         = (*pbTautFlags     & (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) );
    INCHI_MODE bTautFlagsDone     = (*pbTautFlagsDone /*& (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) */);
#if( bRELEASE_VERSION == 0 )
    int bExtract = 0; /*  EXTR_HAS_ATOM_WITH_DEFINED_PARITY; */
#endif

/*^^^ */
    int bFixIsoFixedH = 0;
    int bFixTermHChrg = 0;

#if( TEST_RENUMB_ATOMS == 1 )
    long ulNormTime=0;
    long ulCanonTime=0, ulCanonTime2=0;

    inchiTime ulNormTimeStart;
    inchiTime ulCanonTimeStart;

    InchiTimeGet( &ulNormTimeStart );
#endif

    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    int vABParityUnknown = AB_PARITY_UNDF;
    if ( 0 != ( nUserMode & REQ_MODE_DIFF_UU_STEREO) ) 
    {
        /* Make labels for unknown and undefined stereo different */
        vABParityUnknown = AB_PARITY_UNKN;
    }

    

/*^^^ */
#if( FIX_ISO_FIXEDH_BUG == 1 )
    if (TG_FLAG_FIX_ISO_FIXEDH_BUG & *pbTautFlags)
        bFixIsoFixedH = 1;
#endif
#if( FIX_TERM_H_CHRG_BUG == 1 )
    if (TG_FLAG_FIX_TERM_H_CHRG_BUG & *pbTautFlags)
        bFixTermHChrg = 1;
#endif
/*^^^ */

    memset( s, 0, sizeof(s) );
    if ( pBCN ) {
        memset( pBCN, 0, sizeof( pBCN[0] ) );
    }
    memset( t_group_info, 0, sizeof(*t_group_info) );
    memset( t_group_info_orig, 0, sizeof(*t_group_info_orig) );
    /*norm_at[TAUT_NON] = out_norm_data[TAUT_NON]->at; *//* output normalized non-tautomeric component */
    /*norm_at[TAUT_YES] = out_norm_data[TAUT_YES]->at; *//* output normalized tautomeric component */
    /*norm_at_fixed_bonds[TAUT_NON] = NULL;*/
    /*norm_at_fixed_bonds[TAUT_YES] = out_norm_data[TAUT_YES]->at_fixed_bonds;*/
    for ( i = 0; i < TAUT_NUM; i ++ ) {
        if ( out_norm_data[i]->at ) {
            if ( !(at[i] = (sp_ATOM  *)inchi_malloc( num_inp_at * sizeof(*at[0]) ) ) ) {
                ret = -1;
            }
        } else {
            at[i] = NULL;
        }
    }
    if ( !out_norm_data[TAUT_NON]->at && !out_norm_data[TAUT_YES]->at || !inp_at || ret ) {
        ret = -1;
        goto exit_function;
    }
    /* the first struct to process: tautomeric if exists else non-tautomeric */
    out_at = out_norm_data[TAUT_YES]->at? out_norm_data[TAUT_YES]->at : out_norm_data[TAUT_NON]->at;
    /* copy the input structure to be normalized to the buffer for the normalization data */
    memcpy( out_at, inp_at, num_inp_at*sizeof(out_at[0]) );

    /*  tautomeric groups setting */
    t_group_info->bIgnoreIsotopic = 0;   /*  include tautomeric group isotopic info in MarkTautomerGroups() */
    t_group_info->bTautFlags      = *pbTautFlags;
    t_group_info->bTautFlagsDone  = *pbTautFlagsDone;

    /*  Preprocess the structure; here THE NUMBER OF ATOMS MAY BE REDUCED */
    /*  ??? Ambiguity: H-D may become HD or DH (that is, H+implicit D or D+implicit H) */
    if ( TG_FLAG_H_ALREADY_REMOVED & bTautFlags ) {
        INP_ATOM_DATA *out_norm_data1 = out_norm_data[TAUT_YES]->at? out_norm_data[TAUT_YES] :
                                        out_norm_data[TAUT_NON]->at? out_norm_data[TAUT_NON] : NULL;
        if ( out_norm_data1 ) {
            num_at_tg     =
            num_atoms     = out_norm_data1->num_at - out_norm_data1->num_removed_H;
            num_removed_H = out_norm_data1->num_removed_H;
            t_group_info->tni.nNumRemovedExplicitH = num_removed_H;
        } else {
            ret = -1;
            goto exit_function;
        }
    } else {
        num_at_tg =
        num_atoms = remove_terminal_HDT( num_inp_at, out_at, bFixTermHChrg );
        num_removed_H = num_inp_at - num_atoms;
        t_group_info->tni.nNumRemovedExplicitH = num_removed_H;
        add_DT_to_num_H( num_atoms, out_at );
    }
    /*fix_odd_things( num_atoms, out_at );*/
#if( FIND_RING_SYSTEMS == 1 )
    MarkRingSystemsInp( out_at, num_atoms, 0 );
#endif
    /*  duplicate the preprocessed structure so that all supplied out_norm_data[]->at buffers are filled */
    if ( out_at != out_norm_data[TAUT_YES]->at && out_norm_data[TAUT_YES]->at ) {
        memcpy( out_norm_data[TAUT_YES]->at, out_at, num_inp_at*sizeof(out_at[0]) );
    }
    if ( out_norm_data[TAUT_YES]->at_fixed_bonds && out_norm_data[TAUT_YES]->at ) {
        memcpy( out_norm_data[TAUT_YES]->at_fixed_bonds, out_at, num_inp_at*sizeof(out_at[0]) );
    }
    if ( out_at != out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_NON]->at ) {
        memcpy( out_norm_data[TAUT_NON]->at, out_at, num_inp_at*sizeof(out_at[0]) );
    }

    /*******************************************************************************
     * ??? not true ??? duplicate inp_at and keep inp_at[] unchanged after terminal hydrogens removal
     * set stereo parities in taut_at[], non_taut_at[]
     * obtain max. lenghts of the name stereo parts
     * Ignore absence/presence of isotopic stereo for now
     * mark isotopic atoms
     *******************************************************************************/
    if ( out_norm_data[TAUT_YES]->at && at[TAUT_YES] ) {
        /* final normalization of possibly tautomeric structure */
        ret = mark_alt_bonds_and_taut_groups ( out_norm_data[TAUT_YES]->at, out_norm_data[TAUT_YES]->at_fixed_bonds, num_atoms,
                                               t_group_info, NULL, NULL );
        if ( ret < 0 ) {
            goto exit_function;/*  out of RAM or other normalization problem */
        }
        num_taut_at = ret; /* number of atoms without removed H? */
        num_removed_H_taut = t_group_info->tni.nNumRemovedExplicitH;
        out_norm_data[TAUT_YES]->num_at              = num_atoms + num_removed_H_taut; /* protons might have been removed */
        out_norm_data[TAUT_YES]->num_removed_H       = num_removed_H_taut;
        out_norm_data[TAUT_YES]->nNumRemovedProtons += t_group_info->tni.nNumRemovedProtons;
        for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
            out_norm_data[TAUT_YES]->nNumRemovedProtonsIsotopic[i] += t_group_info->tni.nNumRemovedProtonsIsotopic[i] /*+ t_group_info->num_iso_H[i]*/;
            out_norm_data[TAUT_YES]->num_iso_H[i]                  += t_group_info->num_iso_H[i];
        }
        /* mark deleted isolated tautomeric H(+) */
        if ( num_taut_at == 1 && out_norm_data[TAUT_YES]->at[0].at_type == ATT_PROTON &&
             t_group_info && t_group_info->tni.nNumRemovedProtons == 1 ) {
            out_norm_data[TAUT_YES]->bDeleted = 1;
            FreeInpAtom( &out_norm_data[TAUT_YES]->at_fixed_bonds );
        } else
        if ( (t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) &&
             out_norm_data[TAUT_YES]->at_fixed_bonds) {
             out_norm_data[TAUT_YES]->bTautPreprocessed = 1;
        }
        /*
        if ( !(t_group_info->tni.bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT & ~FLAG_PROTON_SINGLE_REMOVED)) &&
             out_norm_data[TAUT_YES]->at_fixed_bonds) {
             FreeInpAtom( &out_norm_data[TAUT_YES]->at_fixed_bonds );
        }
        */
        /*out_norm_data[TAUT_YES]->num_removed_H = num_removed_H_taut;*/
        out_norm_data[TAUT_YES]->bTautFlags     = *pbTautFlags     = t_group_info->bTautFlags;
        out_norm_data[TAUT_YES]->bTautFlagsDone = *pbTautFlagsDone = t_group_info->bTautFlagsDone;
        out_norm_data[TAUT_YES]->bNormalizationFlags = t_group_info->tni.bNormalizationFlags;
        /* create internal sp_ATOM at[] out of out_norm_data[]->at */
        inp2spATOM( out_norm_data[TAUT_YES]->at, num_inp_at, at[TAUT_YES] );
        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
        ret = set_stereo_parity( out_norm_data[TAUT_YES]->at, at[TAUT_YES], num_taut_at, num_removed_H_taut,
                                 &s[TAUT_YES].nMaxNumStereoAtoms, &s[TAUT_YES].nMaxNumStereoBonds, nUserMode,
                                 bPointedEdgeStereo, vABParityUnknown );
#if( bRELEASE_VERSION == 0 )
        if ( 0 < ret ) {
            bExtract |= EXTR_HAS_ATOM_WITH_DEFINED_PARITY;
        }
        if ( t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT ) {
            bExtract |= EXTR_TAUT_TREATMENT_CHARGES;
        }
#endif
        if ( RETURNED_ERROR( ret ) ) {
            goto exit_function; /*  stereo bond error */
        }
        s[TAUT_YES].bMayHaveStereo    = (s[TAUT_YES].nMaxNumStereoAtoms || s[TAUT_YES].nMaxNumStereoBonds);
        /* 
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */
        s[TAUT_YES].num_isotopic_atoms = set_atom_iso_sort_keys( num_taut_at, at[TAUT_YES], t_group_info,
                                                          &s[TAUT_YES].bHasIsotopicTautGroups );
        /**************************************************************************
         *  prepare tautomeric (if no tautomerism found then prepare non-tautomeric)
         *  structure for canonicalizaton:
         **************************************************************************
         *   remove t-groups that have no H,
         *   remove charges from t-groups if requested
         *   renumber t-groups and find final t_group_info->num_t_groups
         *   add to t-groups lists of endpoints tgroup->nEndpointAtomNumber[]
         *   calculate length of the t-group part of the connection table
         **************************************************************************/
        s[TAUT_YES].nLenLinearCTTautomer = CountTautomerGroups( at[TAUT_YES], num_taut_at, t_group_info );
        if ( RETURNED_ERROR(s[TAUT_YES].nLenLinearCTTautomer) ) { /* added error treatment 9-11-2003 */
            ret = s[TAUT_YES].nLenLinearCTTautomer;
            goto exit_function;
            /*  error has happened; no breakpoint here
            s[TAUT_YES].nLenLinearCTTautomer = 0;
            */
        } else
        if ( s[TAUT_YES].nLenLinearCTTautomer > 0 ) {
            num_at_tg = num_taut_at+t_group_info->num_t_groups;
            /*  ??? -not true- create t_group_info_orig for multiple calls with atom renumbering */
            make_a_copy_of_t_group_info( t_group_info_orig /* dest*/, t_group_info /* source*/ );
            /*  mark isotopic tautomer groups: calculate t_group->iWeight */
            s[TAUT_YES].nLenLinearCTIsotopicTautomer=set_tautomer_iso_sort_keys( t_group_info );
            if ( s[TAUT_YES].nLenLinearCTIsotopicTautomer < 0 ) {
                /* ??? -error cannot happen- error has happened; no breakpoint here */
                s[TAUT_YES].nLenLinearCTIsotopicTautomer = 0;
            }
            out_norm_data[TAUT_YES]->bTautomeric = s[TAUT_YES].nLenLinearCTTautomer;
        }
        /*  new variable: s[TAUT_YES].nLenCT introduced 7-22-2002 */
        GetCanonLengths( num_taut_at, at[TAUT_YES], &s[TAUT_YES], t_group_info );
    }
    if ( out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at && at[TAUT_NON] && !s[TAUT_YES].nLenLinearCTTautomer ) {
        /* the structure is non-tautomeric: use tautomeric treatment results only for it */
        inchi_free( at[TAUT_NON] );
        at[TAUT_NON] = NULL;
    } else
    if ( !out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at &&
         !at[TAUT_NON] && at[TAUT_YES] && !s[TAUT_YES].nLenLinearCTTautomer ) {
        /* requested tautomeric; found non-tautomeric; it is located in out_norm_data[TAUT_YES]->at */
        out_norm_data[TAUT_YES]->bTautomeric = 0;
    } else
    if ( out_norm_data[TAUT_NON]->at && at[TAUT_NON] ) {
        /* the structure needs non-tautomeric treatment: final normalization of non-tautomeric structure */
        ret = mark_alt_bonds_and_taut_groups ( out_norm_data[TAUT_NON]->at, NULL, num_atoms, NULL, &bTautFlags, &bTautFlagsDone );
        if ( ret < 0 ) {
            goto exit_function;  /*  out of RAM or other normalization problem */
        }
        out_norm_data[TAUT_NON]->num_at        = num_atoms + num_removed_H;
        out_norm_data[TAUT_NON]->num_removed_H = num_removed_H;
        out_norm_data[TAUT_NON]->bTautFlags     = *pbTautFlags;
        out_norm_data[TAUT_NON]->bTautFlagsDone = *pbTautFlagsDone;
        out_norm_data[TAUT_NON]->bNormalizationFlags = 0;
        /* create internal sp_ATOM at[] out of out_norm_data[]->at */
        inp2spATOM( out_norm_data[TAUT_NON]->at, num_inp_at, at[TAUT_NON] );
        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
        ret = set_stereo_parity( out_norm_data[TAUT_NON]->at, at[TAUT_NON], num_atoms, num_removed_H,
                                 &s[TAUT_NON].nMaxNumStereoAtoms, &s[TAUT_NON].nMaxNumStereoBonds, nUserMode,
                                 bPointedEdgeStereo, vABParityUnknown );
#if( bRELEASE_VERSION == 0 )
        if ( 0 < ret ) {
            bExtract |= EXTR_HAS_ATOM_WITH_DEFINED_PARITY;
        }
#endif
        if ( RETURNED_ERROR( ret ) ) {
            goto exit_function; /*  stereo bond error */
        }
        s[TAUT_NON].bMayHaveStereo = (s[TAUT_NON].nMaxNumStereoAtoms || s[TAUT_NON].nMaxNumStereoBonds);
        /* 
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */
        s[TAUT_NON].num_isotopic_atoms = set_atom_iso_sort_keys( num_atoms, at[TAUT_NON], NULL, NULL );
        GetCanonLengths( num_atoms, at[TAUT_NON], &s[TAUT_NON], NULL);
        out_norm_data[TAUT_NON]->bTautomeric = 0;
    }

    /**********************************************************/
    /*  common */
    bMayHaveStereo        = s[TAUT_YES].bMayHaveStereo || s[TAUT_NON].bMayHaveStereo;
    bHasIsotopicAtoms     = s[TAUT_NON].num_isotopic_atoms > 0 || s[TAUT_NON].bHasIsotopicTautGroups > 0 ||
                            s[TAUT_YES].num_isotopic_atoms > 0 || s[TAUT_YES].bHasIsotopicTautGroups > 0 ;
/*^^^ */
    if (bFixIsoFixedH)  /* 2008-03-21 DT */
        bHasIsotopicAtoms     = bHasIsotopicAtoms   || 
                                s[TAUT_YES].nLenLinearCTTautomer > 0 && t_group_info &&
                                (0 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[0] ||
                                 1 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[1] ||
                                 2 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[2]) ;
/*^^^ */
    bHasIsotopicAtoms   =   bHasIsotopicAtoms || 
                            s[TAUT_YES].nLenIsotopicEndpoints > 1 && t_group_info &&
                            (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE));

    /*  default mode */
    if ( !(nUserMode & REQ_MODE_DEFAULT) ) {
        /*  default */
        nUserMode |= REQ_MODE_DEFAULT;
    }
    
    /*  adjust the mode to the reality */
    if ( ( nUserMode & REQ_MODE_ISO ) && !bHasIsotopicAtoms ) {
        nUserMode ^= REQ_MODE_ISO;
        nUserMode |= REQ_MODE_NON_ISO;  /*  at least one is needed */
    }
    if ( (nUserMode & REQ_MODE_STEREO) && ( nUserMode & REQ_MODE_ISO ) ) {
        nUserMode |= REQ_MODE_ISO_STEREO;
    }
    if ( (nUserMode & REQ_MODE_STEREO) && !( nUserMode & REQ_MODE_NON_ISO ) ) {
        nUserMode ^= REQ_MODE_STEREO;
    }
    if ( !bMayHaveStereo ) {
        if ( nUserMode & REQ_MODE_STEREO )
            nUserMode ^= REQ_MODE_STEREO;
        if ( nUserMode & REQ_MODE_ISO_STEREO )
            nUserMode ^= REQ_MODE_ISO_STEREO;
    }

    if ( (nUserMode & REQ_MODE_BASIC) && (!out_norm_data[TAUT_NON]->at || !ppINChI[TAUT_NON] || !ppINChI_Aux[TAUT_NON] || !at[TAUT_NON]) ) {
        nUserMode ^= REQ_MODE_BASIC;
    }
    if ( (nUserMode & REQ_MODE_TAUT) && (!out_norm_data[TAUT_YES]->at || !ppINChI[TAUT_YES] || !ppINChI_Aux[TAUT_YES] || !at[TAUT_YES]) ) {
        nUserMode ^= REQ_MODE_TAUT;
    }
        

    switch ((int)nUserMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) {
    case REQ_MODE_BASIC:
        n1 = TAUT_NON;
        n2 = TAUT_NON;
        break;
    case REQ_MODE_TAUT:
        n1 = TAUT_YES;
        n2 = TAUT_YES;
        break;
    case (REQ_MODE_BASIC | REQ_MODE_TAUT):
        n1 = TAUT_NON;
        n2 = TAUT_YES;
        break;
    default:
        ret = -3;
        goto exit_function; /*  program error: inconsistent nUserMode or missing taut/non-taut allocation */ /*   <BRKPT> */
    }
#if( TEST_RENUMB_ATOMS == 1 )
    ulNormTime = InchiTimeElapsed( &ulNormTimeStart);
#endif
    /************************************************************
     *                                                          *
     *       Obtain all non-stereo canonical numberings         *
     *                                                          *
     ************************************************************/
#if( TEST_RENUMB_ATOMS == 1 )
        InchiTimeGet( &ulCanonTimeStart );
#endif
    if ( (nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO) ) {
        /* added for special non-isotopic test mode 2004-10-04 */
        if ( t_group_info ) {
            t_group_info->bIgnoreIsotopic = 1;
            if ( t_group_info->nIsotopicEndpointAtomNumber ) {
                t_group_info->nIsotopicEndpointAtomNumber[0] = inchi_min(1, t_group_info->nIsotopicEndpointAtomNumber[0]);
            }
            memset( t_group_info->num_iso_H, 0, sizeof(t_group_info->num_iso_H) );
            memset ( t_group_info->tni.nNumRemovedProtonsIsotopic, 0, sizeof(t_group_info->tni.nNumRemovedProtonsIsotopic));
            t_group_info->bTautFlagsDone &= ~(TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
        }
        for ( i = 0; i < TAUT_NUM; i ++ ) {
            s[i].bHasIsotopicTautGroups = 0;
            s[i].bIgnoreIsotopic = 1;
            s[i].nLenIsotopic = 0;
            s[i].nLenIsotopicEndpoints = 0;
            s[i].nLenLinearCTIsotopicTautomer = 0;
            s[i].num_isotopic_atoms = 0;
        }
        bHasIsotopicAtoms = 0;
    }
    ret = GetBaseCanonRanking( num_atoms, num_at_tg, at, t_group_info, s, pBCN, ulMaxTime, bFixIsoFixedH );
#if( TEST_RENUMB_ATOMS == 1 )
        ulCanonTime = InchiTimeElapsed( &ulCanonTimeStart );
#endif
    if ( ret < 0 ) {
        goto exit_function; /*  program error */
    }
#if( bRELEASE_VERSION == 0 && FIND_CANON_NE_EQUITABLE == 1 )
    /* Debug only: find whether canonical equivalence is different from equitable partition */
    if ( bCanonIsFinerThanEquitablePartition( num_atoms, at[n1], pBCN->ftcn[TAUT_NON].nSymmRankCt ) ) {
        bExtract |= EXTR_CANON_NE_EQUITABLE;
    }
#endif
    /* added for special non-isotopic test mode 2004-10-04 */
    if ( !pBCN->ftcn[n1].PartitionCt.Rank ) {
        n1 = ALT_TAUT(n1);
    }
    if ( !pBCN->ftcn[n2].PartitionCt.Rank ) {
        n2 = ALT_TAUT(n2);
    }
    if ( n1 > n2 ) {
        ret = CT_TAUCOUNT_ERR;
        goto exit_function; /*  program error */
    }

    /************************************************************
     *                                                          *
     *       Obtain stereo canonical numberings                 *
     *                                                          *
     ************************************************************/

    for ( i = n2; i >= n1 && !RETURNED_ERROR( ret ); i -- ) {

        memset( pCS, 0, sizeof(*pCS) );

        switch( i ) {
        case TAUT_NON: /*  non-tautomeric */
            nMode  = 0;
            nMode  = (s[i].nLenLinearCTTautomer == 0)? CANON_MODE_CT:CANON_MODE_TAUT;
            nMode |= (bHasIsotopicAtoms && (nUserMode & REQ_MODE_ISO))? CANON_MODE_ISO:0;
            nMode |= (s[TAUT_NON].bMayHaveStereo && (nUserMode & REQ_MODE_STEREO) )? CANON_MODE_STEREO:0;
            nMode |= (bHasIsotopicAtoms && s[TAUT_NON].bMayHaveStereo && (nUserMode & REQ_MODE_ISO_STEREO))? CANON_MODE_ISO_STEREO:0;
            nMode |= (nUserMode & REQ_MODE_NOEQ_STEREO   )? CMODE_NOEQ_STEREO    : 0;
            nMode |= (nUserMode & REQ_MODE_REDNDNT_STEREO)? CMODE_REDNDNT_STEREO : 0;
            nMode |= (nUserMode & REQ_MODE_NO_ALT_SBONDS )? CMODE_NO_ALT_SBONDS  : 0;
            if ( (nMode & CANON_MODE_STEREO)     == CANON_MODE_STEREO ||
                 (nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO ) {
                nMode |= (nUserMode & REQ_MODE_RELATIVE_STEREO)? CMODE_RELATIVE_STEREO: 0;
                nMode |= (nUserMode & REQ_MODE_RACEMIC_STEREO )? CMODE_RACEMIC_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_SC_IGN_ALL_UU  )? CMODE_SC_IGN_ALL_UU  : 0;
                nMode |= (nUserMode & REQ_MODE_SB_IGN_ALL_UU  )? CMODE_SB_IGN_ALL_UU  : 0;
            }
            if ( ret= AllocateCS( pCS, num_atoms, num_atoms, s[TAUT_NON].nLenCT, s[TAUT_NON].nLenCTAtOnly,
                             s[TAUT_NON].nLenLinearCTStereoDble, s[TAUT_NON].nMaxNumStereoBonds,
                             s[TAUT_NON].nLenLinearCTStereoCarb, s[TAUT_NON].nMaxNumStereoAtoms,
                             0, 0, s[TAUT_NON].nLenIsotopic, nMode, pBCN ) ) {
                goto exit_function;
            }
            *pCS2 = *pCS;
            break;
        case TAUT_YES: /*  tautomeric */
            nMode  = 0;
            nMode  = (s[i].nLenLinearCTTautomer == 0)? CANON_MODE_CT:CANON_MODE_TAUT;
            nMode |= (bHasIsotopicAtoms && (nUserMode & REQ_MODE_ISO) )? CANON_MODE_ISO:0;
            nMode |= (s[TAUT_YES].bMayHaveStereo && (nUserMode & REQ_MODE_STEREO) )? CANON_MODE_STEREO:0;
            nMode |= (bHasIsotopicAtoms && s[TAUT_YES].bMayHaveStereo && (nUserMode & REQ_MODE_ISO_STEREO))? CANON_MODE_ISO_STEREO:0;
            nMode |= (nUserMode & REQ_MODE_NOEQ_STEREO   )? CMODE_NOEQ_STEREO    : 0;
            nMode |= (nUserMode & REQ_MODE_REDNDNT_STEREO)? CMODE_REDNDNT_STEREO : 0;
            nMode |= (nUserMode & REQ_MODE_NO_ALT_SBONDS )? CMODE_NO_ALT_SBONDS  : 0;
            if ( (nMode & CANON_MODE_STEREO)     == CANON_MODE_STEREO ||
                 (nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO ) {
                nMode |= (nUserMode & REQ_MODE_RELATIVE_STEREO)? CMODE_RELATIVE_STEREO: 0;
                nMode |= (nUserMode & REQ_MODE_RACEMIC_STEREO )? CMODE_RACEMIC_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_SC_IGN_ALL_UU  )? CMODE_SC_IGN_ALL_UU  : 0;
                nMode |= (nUserMode & REQ_MODE_SB_IGN_ALL_UU  )? CMODE_SB_IGN_ALL_UU  : 0;
            }
            if ( ret= AllocateCS( pCS, num_atoms, num_at_tg, s[TAUT_YES].nLenCT, s[TAUT_YES].nLenCTAtOnly,
                             s[TAUT_YES].nLenLinearCTStereoDble, s[TAUT_YES].nMaxNumStereoBonds,
                             s[TAUT_YES].nLenLinearCTStereoCarb, s[TAUT_YES].nMaxNumStereoAtoms,
                             s[TAUT_YES].nLenLinearCTTautomer, s[TAUT_YES].nLenLinearCTIsotopicTautomer,
                             s[TAUT_YES].nLenIsotopic, nMode, pBCN ) ) {
                goto exit_function;
            }
            *pCS2 = *pCS;
            break;
        }

    
        /*^^^ 2009-12-05 */
        nMode |= (nUserMode & REQ_MODE_DIFF_UU_STEREO)? REQ_MODE_DIFF_UU_STEREO : 0;
        /*^^^ 2009-12-05 */


        /*  settings */
        pCS->lNumDecreasedCT               = -1;                             
        pCS->bDoubleBondSquare             =  DOUBLE_BOND_NEIGH_LIST? 2:0;  /*  2 => special mode */
        pCS->bIgnoreIsotopic               =  !((s[TAUT_NON].num_isotopic_atoms ||        
                                                s[TAUT_YES].num_isotopic_atoms ||    
                                                s[TAUT_YES].bHasIsotopicTautGroups) ||
                                               (nUserMode & REQ_MODE_NON_ISO) ||
                                               !(nUserMode & REQ_MODE_ISO));
        
        if ( (nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO) ) {
            pCS->bIgnoreIsotopic = 1; /* 10-04-2004 */
        }
        
        if ( i == TAUT_YES ) {  /* tautomeric */
            pCS->t_group_info                  = t_group_info; /*  ??? make a copy or reuse ???  */
            pCS->t_group_info->bIgnoreIsotopic = !(s[TAUT_YES].bHasIsotopicTautGroups ||
                                                   (nUserMode & REQ_MODE_NON_ISO) ||
                                                   !(nUserMode & REQ_MODE_ISO));
            if ( (nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO) ) {
                pCS->t_group_info->bIgnoreIsotopic = 1; /* 10-04-2004 */
            }
        }
        pCS->ulTimeOutTime  = pBCN->ulTimeOutTime;
        /*=========== Obsolete Mode Bits (bit 0 is Least Significant Bit) ===========
         *
         *  Mode      Bits       Description                                
         *   '0' c    0          Only one connection table canonicalization 
         *   '1' C    1          Recalculate CT using fixed nSymmRank       
         *   '2' i    1|2        Isotopic canonicalization (internal)       
         *   '3' I    1|2|4      Isotopic canonicalization (output)
         *   '4' s    1|8        Stereo canonicalization                    
         *   '5' S    1|2|4|16   Stereo isotopic canonicalization           
         *   '6' A    1|2|4|8|16 Output All
         */
#if( TEST_RENUMB_ATOMS == 1 )
        InchiTimeGet( &ulCanonTimeStart );
#endif
        /***************************************
           The last canonicalization step
         ***************************************/
        if ( pBCN ) {
            /* USE_CANON2 == 1 */
            pCS->NeighList  = NULL;
            pCS->pBCN       = pBCN;
            ret = Canon_INChI( num_atoms, i?num_at_tg:num_atoms, at[i], pCS, nMode, i);
        } else {
            /* old way */
            pCS->NeighList  = CreateNeighList( num_atoms, i?num_at_tg:num_atoms, at[i], pCS->bDoubleBondSquare, pCS->t_group_info );
            pCS->pBCN       = NULL;
            ret = Canon_INChI( num_atoms, i?num_at_tg:num_atoms, at[i], pCS, nMode, i);
        }

        pINChI     = ppINChI[i];      /* pointers to already allocated still empty InChI */
        pINChI_Aux = ppINChI_Aux[i];
        if ( ret <= 0 ) {
            /***************************************/
            /*  failure in Canon_INChI()            */
            /***************************************/
            pINChI->nErrorCode     = ret;
            pINChI_Aux->nErrorCode = ret;
        } else {
            /***************************************/
            /*  success Canon_INChI()               */
            /*  save canonicalization results in   */
            /*  pINChI and pINChI_Aux                */
            /***************************************/
            pINChI->nErrorCode     = 0;
            pINChI_Aux->nErrorCode = 0;
            pINChI->bDeleted = pINChI_Aux->bDeleted = out_norm_data[i]->bDeleted;
            pINChI_Aux->nCanonFlags         = pCS->nCanonFlags;
            pINChI_Aux->bTautFlags          = out_norm_data[i]->bTautFlags;
            pINChI_Aux->bTautFlagsDone      = out_norm_data[i]->bTautFlagsDone;
            pINChI_Aux->bNormalizationFlags = out_norm_data[i]->bNormalizationFlags;
            /*  may return an error or a warning */
            ret = FillOutINChI( pINChI, pINChI_Aux,
                               num_atoms, i?num_at_tg:num_atoms,
                               i?num_removed_H_taut:num_removed_H, at[i], 
                               out_norm_data[i]->at, pCS, i, nUserMode, 
                               pStrErrStruct );
            if ( RETURNED_ERROR( ret ) ) {
                /* failure in FillOutINChI() */
                pINChI->nErrorCode      = ret;
                pINChI_Aux->nErrorCode  = ret;
            } else {
                /****************************/
                /* success in FillOutINChI() */
                /****************************/
#if( bRELEASE_VERSION == 0 )
                if ( pINChI->Stereo &&
                     (pINChI->Stereo->nCompInv2Abs && !pINChI->Stereo->bTrivialInv) ||
                     pINChI->StereoIsotopic &&
                     (pINChI->StereoIsotopic->nCompInv2Abs && !pINChI->StereoIsotopic->bTrivialInv) ) {
                    bExtract |= EXTR_NON_TRIVIAL_STEREO;
                }
#endif
                /* mark non-tautomeric representation as having another, tautomeric representation */
                if ( pINChI_Aux && s[TAUT_YES].nLenLinearCTTautomer ) {
                    pINChI_Aux->bIsTautomeric = s[TAUT_YES].nLenLinearCTTautomer;
                }
#if( bRELEASE_VERSION == 0 )
                pCS->bExtract   |= bExtract;
                pINChI->bExtract |= pCS->bExtract;
#endif
                ret2 = CheckCanonNumberingCorrectness(
                                num_atoms, i?num_at_tg:num_atoms,
                                at[i], pCS, i, pStrErrStruct );
                if ( ret2 ) {
                    pINChI->nErrorCode      = ret2;
                    pINChI_Aux->nErrorCode  = ret2;
                    ret = ret2;
                }
            }
        }
#if( TEST_RENUMB_ATOMS == 1 )
        ulCanonTime2 = InchiTimeElapsed( &ulCanonTimeStart );
        pINChI_Aux->ulCanonTime = ulCanonTime+ulCanonTime2;
        pINChI_Aux->ulNormTime  = ulNormTime;
#endif
        FreeNeighList( pCS->NeighList );
        DeAllocateCS( pCS2 );

        pINChI = NULL;      /* avoid dangling pointers */
        pINChI_Aux = NULL;  /* avoid dangling pointers */
    }
    if ( ret == 0 ) {
        ret = num_atoms;
    }
    /*  treat the results later */

exit_function:
    DeAllocBCN( pBCN );
    if ( at[TAUT_YES] )
        inchi_free( at[TAUT_YES] );
    if ( at[TAUT_NON] )
        inchi_free( at[TAUT_NON] );
    if ( ti_out ) {
        *ti_out = *t_group_info;
    } else {
        free_t_group_info( t_group_info );
    }
    free_t_group_info( t_group_info_orig );
    return ret;
}
#ifndef INCHI_ANSI_ONLY /* { */
/***************************************************************************************/
int GetAtomOrdNbrInCanonOrd( inp_ATOM *norm_at, AT_NUMB *nAtomOrdNbr,
                            AT_NUMB *nOrigAtNosInCanonOrd, int num_at )
{
    AT_NUMB   *nCanonNbr, *nOrigAtNos, *nOrigAtNosOrd;
    int          i, ret;

    ret = 0;

    nCanonNbr     = (AT_NUMB *)inchi_calloc( num_at, sizeof(nCanonNbr[0]) );
    nOrigAtNos    = (AT_NUMB *)inchi_calloc( num_at, sizeof(nOrigAtNos[0]) );
    nOrigAtNosOrd = (AT_NUMB *)inchi_calloc( num_at, sizeof(nOrigAtNosOrd[0]) );

    if ( !nCanonNbr || !nOrigAtNos || !nAtomOrdNbr || !nOrigAtNosOrd ) {
        ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
        goto exit_function;
    }
    for ( i = 0; i < num_at; i ++ ) {
        nCanonNbr[i] = nAtomOrdNbr[i] = nOrigAtNosOrd[i] = (AT_NUMB)i;
        nOrigAtNos[i] = norm_at[i].orig_at_number;
    }
    /*  get nCanonNbr[]: canon. numbers-1 in order of increasing original atom numbers */
    pn_RankForSort = nOrigAtNosInCanonOrd;
    qsort( nCanonNbr, num_at, sizeof(nCanonNbr[0]), CompRank );
    /*  get nOrigAtNosOrd[]: norm_atom ord. numbers the same order of increasing original atom numbers */
    pn_RankForSort = nOrigAtNos;
    qsort( nOrigAtNosOrd, num_at, sizeof(nOrigAtNosOrd[0]), CompRank ); 
    /*  check whether the 2 sets of origiginal atom numbers have identical elements */
    for ( i = 0; i < num_at; i ++ ) {
        if ( nOrigAtNosInCanonOrd[nCanonNbr[i]] != nOrigAtNos[nOrigAtNosOrd[i]] ) {
            ret = CT_RANKING_ERR;  /*   <BRKPT> */
            goto exit_function;
        }
    }
    for ( i = 0; i < num_at; i ++ ) {
        nAtomOrdNbr[(int)nCanonNbr[i]] = nOrigAtNosOrd[i];
    }

/*
    pn_RankForSort = nCanonNbr;
    qsort( nAtomOrdNbr, num_at, sizeof(nCanonNbr[0]), CompRank );
*/

exit_function:
    if ( nCanonNbr )
        inchi_free( nCanonNbr );
    if ( nOrigAtNos )
        inchi_free( nOrigAtNos );
    if ( nOrigAtNosOrd )
        inchi_free( nOrigAtNosOrd );
    return ret;
}
/***************************************************************************************/
int FillOutCanonInfAtom(inp_ATOM *norm_at, INF_ATOM_DATA *inf_norm_at_data, int init_num_at, int bIsotopic,
                        INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode)
{
    int          i, j, m, n, num_stereo, k, c, ret, len_str, len, atw, num_err;
    int          next_atom[MAX_CUMULENE_LEN+1], best_next_atom[MAX_CUMULENE_LEN+1], cur_atom;
    int          next_neigh[MAX_CUMULENE_LEN+1], best_next_neigh[MAX_CUMULENE_LEN+1], best_len;
    int          num_iso_H[NUM_H_ISOTOPES];
    char        *str;
    AT_NUMB      g, e;
    int          num_at           = pINChI->nNumberOfAtoms;
    int          nNumberOfTGroups = (pINChI->lenTautomer && pINChI->nTautomer && pINChI->nTautomer[0])? (int)pINChI->nTautomer[0] : 0;
    AT_NUMB     *nOrigAtNosInCanonOrd;
    INChI_Stereo *Stereo;
    AT_NUMB     *nConstitEquNumbers;
    AT_NUMB     *nConstitEquTGroupNumbers;
    S_CHAR      *t_parity = NULL;
    AT_NUMB     *nNumber = NULL;
    int          bIncludeIsotopicH;

    AT_NUMB     *nNormAtNosInCanonOrd;
    int (*MakeNumber)(char*, int, const char*, int) = bAbcNumbers? MakeAbcNumber:MakeDecNumber;
    int bRel= (0 != (nMode & ( REQ_MODE_RELATIVE_STEREO)));
    int bRac= (0 != (nMode & ( REQ_MODE_RACEMIC_STEREO)));
    int bRelRac= bRel || bRac;
    int bDoDisplaySp3 = 1;

    inf_ATOM *inf_norm_at = inf_norm_at_data? inf_norm_at_data->at : NULL;

    ret = 0;
    num_err = 0;

    if ( !inf_norm_at )
        return ret;
    /* prepare removeable protons and H info */
    inf_norm_at_data->nNumRemovedProtons = pINChI_Aux->nNumRemovedProtons;
    MakeRemovedProtonsString( pINChI_Aux->nNumRemovedProtons, pINChI_Aux->nNumRemovedIsotopicH, NULL, bIsotopic,
                              inf_norm_at_data->szRemovedProtons, &inf_norm_at_data->num_removed_iso_H );
    /* fill out info atom */
    if ( bIsotopic && !(pINChI->nNumberOfIsotopicAtoms || pINChI->nNumberOfIsotopicTGroups ||
                        pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0]>1) )
        bIsotopic = 0;
    
    Stereo                   = bIsotopic? pINChI->StereoIsotopic :
                                          pINChI->Stereo;
    bDoDisplaySp3 = (NULL != Stereo) && (Stereo->nNumberOfStereoCenters > 0); 

#if( REL_RAC_STEREO_IGN_1_SC == 1 )
    if ( bDoDisplaySp3 && bRelRac && Stereo->nNumberOfStereoCenters < 2 &&
         (Stereo->nCompInv2Abs || ATOM_PARITY_ILL_DEF(Stereo->t_parity[0]) ) ) {
        bDoDisplaySp3 = 0;
        if ( Stereo->nCompInv2Abs ) {
            inf_norm_at_data->StereoFlags |= bRel? INF_STEREO_REL : bRac? INF_STEREO_RAC : 0;
        }
    }
#endif    
    /* flag has stereo */
    if ( (NULL != Stereo) && (bDoDisplaySp3 || Stereo->nNumberOfStereoBonds > 0) ) {
        inf_norm_at_data->StereoFlags |= INF_STEREO;
    }

    /*
    if ( bDoDisplaySp3 && bRelRac && Stereo->nNumberOfStereoCenters < 2 &&
         (Stereo->nCompInv2Abs || ATOM_PARITY_ILL_DEF(Stereo->t_parity[0]) ) ) {
        bDoDisplaySp3 = 0;
    }
    */
    if ( bDoDisplaySp3 && Stereo->nCompInv2Abs ) {
        /* inversion changes stereo */
        if ( bRel ) {
            inf_norm_at_data->StereoFlags |= INF_STEREO_REL;
        } else
        if ( bRac ) {
            inf_norm_at_data->StereoFlags |= INF_STEREO_RAC;
        } else {
            inf_norm_at_data->StereoFlags |= INF_STEREO_ABS;
        }
        if ( bRelRac ) {
            inf_norm_at_data->StereoFlags |= (Stereo->nCompInv2Abs > 0)? INF_STEREO_NORM : INF_STEREO_INV;
        }
    }
    if ( bDoDisplaySp3 && Stereo->nCompInv2Abs < 0 && !bRelRac ) {
       /* display Inv stereo which is Absolute Stereo */
        nNumber  = Stereo->nNumberInv;
        t_parity = Stereo->t_parityInv;
        nOrigAtNosInCanonOrd     = bIsotopic? pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv :
                                              pINChI_Aux->nOrigAtNosInCanonOrdInv;
    } else {
       /* display Inv stereo which is Absolute Stereo */
        if ( bDoDisplaySp3 ) {
            nNumber  = Stereo->nNumber;
            t_parity = Stereo->t_parity;
        }
        nOrigAtNosInCanonOrd     = bIsotopic? pINChI_Aux->nIsotopicOrigAtNosInCanonOrd :
                                              pINChI_Aux->nOrigAtNosInCanonOrd;
    }

    nConstitEquNumbers       = bIsotopic? pINChI_Aux->nConstitEquIsotopicNumbers :
                                          pINChI_Aux->nConstitEquNumbers;
    nConstitEquTGroupNumbers = bIsotopic? pINChI_Aux->nConstitEquIsotopicTGroupNumbers :
                                          pINChI_Aux->nConstitEquTGroupNumbers;
    memset( inf_norm_at, 0, init_num_at*sizeof(inf_norm_at[0]) );

    /*  obtain norm_at[] atom numbers (from zero) in order of canonical numbers */
    nNormAtNosInCanonOrd = (AT_NUMB *)inchi_calloc( num_at, sizeof(nNormAtNosInCanonOrd[0]) );
    if ( ret = GetAtomOrdNbrInCanonOrd( norm_at, nNormAtNosInCanonOrd, nOrigAtNosInCanonOrd, num_at ) ) {
        goto exit_function;
    }
    
    /*  atom canonical and equivalence numbers > 0 */
    for ( i = 0; i < num_at; i ++ ) {
        j = (int)nNormAtNosInCanonOrd[i];
        if ( j < 0 || j >= num_at )
            continue;
        inf_norm_at[j].nCanonNbr = (AT_NUMB)(i+1);
        inf_norm_at[j].nCanonEquNbr = nConstitEquNumbers[i];
#ifdef DISPLAY_DEBUG_DATA
        inf_norm_at[j].nDebugData = 0;
#if( DISPLAY_DEBUG_DATA == DISPLAY_DEBUG_DATA_C_POINT )
        inf_norm_at[j].nDebugData = norm_at[j].c_point;
#endif
#endif
    }
    /*  tautomeric groups */
    if ( nNumberOfTGroups ) {
        /*
         :   start from 1: bypass number of t-groups
         :   j is a counter within the current t-group
         :   g is a tautomeric group canonical number
         :   e is a tautomeric group equivalence
         */
        for ( g = 1, i = 1; g <= nNumberOfTGroups; g ++ ) {
            n = (int)pINChI->nTautomer[i] - INCHI_T_NUM_MOVABLE; /*  number of atoms in t-group */
            e = nConstitEquTGroupNumbers[(int)g - 1];
            /*  bypass number of hydrogen atoms, negative charges, ... */
            for ( i += INCHI_T_NUM_MOVABLE+1, j = 0; j < n && i < pINChI->lenTautomer; j ++, i ++ ) {
                /*  scan canonical numbers of atoms within the atom t-group */
                k = (int)nNormAtNosInCanonOrd[(int)pINChI->nTautomer[i]-1];
                inf_norm_at[k].nTautGroupCanonNbr = g;
                inf_norm_at[k].nTautGroupEquNbr   = e;
            }
        }
        if ( i != pINChI->lenTautomer || g != nNumberOfTGroups+1 ) {
            ret = CT_TAUCOUNT_ERR;  /*   <BRKPT> */
            goto exit_function;
        }
    }
    /* atoms that may exchange isotopic H */
    if ( bIsotopic && pINChI->nPossibleLocationsOfIsotopicH && (n = (int)pINChI->nPossibleLocationsOfIsotopicH[0]) ) {
        for ( i = 1; i < n; i ++ ) {
            j = (int)pINChI->nPossibleLocationsOfIsotopicH[i];
            k = (int)nNormAtNosInCanonOrd[j - 1];
            if ( !inf_norm_at[k].nTautGroupCanonNbr ) {
                inf_norm_at[k].cFlags |= AT_FLAG_ISO_H_POINT;
            }
        }
    }

#if( DISPLAY_RING_SYSTEMS == 1 )
    /*  debug only */
    for ( j = 0; j < num_at; j ++ ) {
        inf_norm_at[j].nCanonNbr          = norm_at[j].nBlockSystem;
        inf_norm_at[j].nCanonEquNbr       = norm_at[j].nRingSystem;
#if( USE_DISTANCES_FOR_RANKING == 1 )
        inf_norm_at[j].nTautGroupCanonNbr   = norm_at[j].nDistanceFromTerminal;
        inf_norm_at[j].nTautGroupEquNbr     = norm_at[j].bCutVertex;
#else
        inf_norm_at[j].nTautGroupCanonNbr   = norm_at[j].bCutVertex;
        inf_norm_at[j].nTautGroupEquNbr   = 0;
#endif
    }
#endif



    /*  Write isotopic mass, chemical element symbols and hydrogens, charge, radical, canon. numbers */
    len_str = sizeof(inf_norm_at[0].at_string);
    for ( i = 0; i < init_num_at; i ++ ) {
        str = inf_norm_at[i].at_string;
        len = 0;
        bIncludeIsotopicH = bIsotopic && !inf_norm_at[i].nTautGroupCanonNbr && !(inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT);
        /*  isotopic mass */
        atw = 0;
        if ( norm_at[i].iso_atw_diff && bIsotopic ) {
            if ( norm_at[i].at_type == ATT_PROTON ) {
                ; /* do not set isotopic mass of a tautomeric proton */
            } else
            if ( norm_at[i].el_number == PERIODIC_NUMBER_H && norm_at[i].chem_bonds_valence == 1 &&
                 !norm_at[i].charge && !norm_at[i].radical && !norm_at[i].num_H &&
                 (inf_norm_at[j=(int)norm_at[i].neighbor[0]].nTautGroupCanonNbr || (inf_norm_at[j].cFlags & AT_FLAG_ISO_H_POINT) ) ) {
                ; /* do not set isotopic mass of an exchangeable proton */
            } else {
                atw = get_atw(norm_at[i].elname);
                atw += (norm_at[i].iso_atw_diff>0)? norm_at[i].iso_atw_diff-1:norm_at[i].iso_atw_diff;
                /*len += sprintf( str+len, "^%d", atw );*/
            }
        }
        /*  element name */
        if ( norm_at[i].el_number == PERIODIC_NUMBER_H && 2 <= atw && atw <= 3 ) {
            len += sprintf( str+len, "%s", atw==2? "D" : "T" );
        } else {
            if ( atw ) {
                len += sprintf( str+len, "^%d", atw );
            }
            len += sprintf( str+len, "%s", norm_at[i].elname );
        }
        /*  hydrogens */
        /*  find number of previuosly removed terminal hydrogen atoms because these terminal H will be displayed */
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
            num_iso_H[j] = norm_at[i].num_iso_H[j];
        }
        /* n = number of implicit H to display */
        for ( j = num_at, n = (int)norm_at[i].num_H; j < init_num_at; j ++ ) {
            /*  subtract number of removed terminal */
            /*  H atoms from the total number of H atoms */
            if ( i == (int)norm_at[j].neighbor[0] ) {
                n -= 1; /* found explicit H => decrement number of implicit H */
                m = (int)norm_at[j].iso_atw_diff-1;
                if ( 0 <= m && m < NUM_H_ISOTOPES ) {
                    /*  subtract number of removed terminal isotopic H */
                    /*  atoms from the total number of isotopic H atoms */
                    num_iso_H[m] -= 1;
                }
            }
        }
        /* at this point n = number of implicit H to display, 
           num_iso_H[] contains number of implicit isotopic H among n */
        if ( bIncludeIsotopicH ) {
            /*  subtract number of isotopic H atoms from the total number of H atoms */
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                n -= num_iso_H[j];
            }
        }
        /*  non-isotopic hydrogen atoms */
        if ( n > 1 ) {
            len += sprintf( str+len, "H%d", n );
        } else
        if ( n == 1 ) {
            len += sprintf( str+len, "H" );
        }
        /*  isotopic hydrogen atoms */
        if ( bIncludeIsotopicH ) {
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                if ( num_iso_H[j] ) {
                    if ( j == 0 || j != 1 && j != 2 ) {
                        len += sprintf( str+len, "^%dH", j+1 );
                    } else {
                        len += sprintf( str+len, j == 1? "D" : "T" );
                    }
                    if ( num_iso_H[j] != 1 ) {
                        len += sprintf( str+len, "%d", (int)num_iso_H[j] );
                    }
                }
            }
        }
        if ( norm_at[i].el_number == PERIODIC_NUMBER_H && str[0] == str[1] ) {
            char *q;
            if ( !str[2] ) {
                str[1] = '2';  /* quick fix: replace HH with H2 */
            } else
            if ( isdigit( UCINT str[2] ) && (n = strtol( str+2, &q, 10 )) && !q[0] ) {
                len = 1 + sprintf( str+1, "%d", n+1 );
            }
        }
        /*
        if (  str[0] == 'H' && str[1] == 'H' && !str[2] ) {
            str[1] = '2';
        }
        */
        /*  charge */
        if ( abs(norm_at[i].charge) > 1 )
            len += sprintf( str+len, "%+d", norm_at[i].charge );
        else
        if ( abs(norm_at[i].charge) == 1 )
            len += sprintf( str+len, "%s", norm_at[i].charge>0? "+" : "-" );
        /*  radical */
        if ( norm_at[i].radical )
            len += sprintf( str+len, "%s", norm_at[i].radical==RADICAL_SINGLET? ":" :
                                           norm_at[i].radical==RADICAL_DOUBLET? "." :
                                           norm_at[i].radical==RADICAL_TRIPLET? ".." : "?");
    }

    /*  stereogenic centers */
    if ( bDoDisplaySp3 && Stereo && 0 < (num_stereo = Stereo->nNumberOfStereoCenters) ) {
        for ( i = 0; i < num_stereo; i ++ ) {
            j = (int)nNormAtNosInCanonOrd[(int)nNumber[i]-1];
            c = t_parity[i];
            c = c==1? '-' : c==2? '+' : c==3? 'u' : c== 4? '?' : '*';
            inf_norm_at[j].cStereoCenterParity = c;
            str=inf_norm_at[j].at_string;
            len = strlen(str);
            if ( len + 3 < (int)sizeof(inf_norm_at[0].at_string) ) {
                str[len++] = '(';
                str[len++] = inf_norm_at[j].cStereoCenterParity;
                str[len++] = ')';
                str[len] = '\0';
                /*  mark ambuguous stereo center */
                if ( norm_at[j].bAmbiguousStereo && (c=='+' || c=='-' || c=='?') &&  str[0] != '!' &&
                     len+1 < (int)sizeof(inf_norm_at[0].at_string) ) {
                    memmove( str+1, str, len+1 );
                    str[0] = '!'; /* output the atom in red color */
                }
            }
        }
    }
    
    /*  stereogenic bonds */
    /*  (cumulenes with odd number of double bonds are stereocenters, */
    /*   and atom parity should be set) */
    if ( Stereo && 0 < (num_stereo = Stereo->nNumberOfStereoBonds) ) {
        for ( i = 0; i < num_stereo; i ++ ) {
            int start_at, num_eql=0, bAmbiguousStereoBond = 0;
            j = (int)nNormAtNosInCanonOrd[(int)Stereo->nBondAtom1[i]-1];
            k = (int)nNormAtNosInCanonOrd[(int)Stereo->nBondAtom2[i]-1];
            start_at = j;
            c = Stereo->b_parity[i];
                                
            c = c==1? '-' : c==2? '+' : c==3? 'u' : c== 4? '?' : '*';

            /*  mark ambuguous stereo bond atom(s) */
            if ( norm_at[j].bAmbiguousStereo && (c=='+' || c=='-' ) &&
                 (len=strlen(str=inf_norm_at[j].at_string)+1) < (int)sizeof(inf_norm_at[0].at_string) &&
                  str[0] != '!' ) {
                memmove( str+1, str, len );
                str[0] = '!'; /* output the atom in red color */
                bAmbiguousStereoBond ++;
            }
            if ( norm_at[k].bAmbiguousStereo && (c=='+' || c=='-') &&
                 (len=strlen(str=inf_norm_at[k].at_string)+1) < (int)sizeof(inf_norm_at[0].at_string) &&
                  str[0] != '!' ) {
                memmove( str+1, str, len );
                str[0] = '!'; /* output the atom in red color */
                bAmbiguousStereoBond ++;
            }

            /*  find the opposite atom k. */
            /*  Note: since it may be a cumulene, find the shortest(best) path */
            /*  to atom number k to avoid confusion in case of, for example, */
            /*  4-member aromatic rings. */
            best_len = MAX_CUMULENE_LEN+1; /* moved here from inside the cycle 1-8-2003 */
            for ( n = 0; n < norm_at[j].valence; n ++ ) {
                if ( norm_at[j].bond_type[n] == BOND_SINGLE ) {
                    /*  single bond cannot be stereogenic. */
                    continue;
                }
                /* best_len = MAX_CUMULENE_LEN+1; */
                len = 0; /*  number of bonds in cumulene - 1 */
                next_atom[len]  = (int)norm_at[j].neighbor[n];
                next_neigh[len] = n;
                cur_atom = j;
                while ( next_atom[len] != k && len < MAX_CUMULENE_LEN && 2 == norm_at[next_atom[len]].valence ) {
                    next_neigh[len+1] = ((int)norm_at[next_atom[len]].neighbor[0] == cur_atom);
                    next_atom[len+1]  =  (int)norm_at[next_atom[len]].neighbor[next_neigh[len+1]];
                    cur_atom = next_atom[len];
                    len ++;
                }
                if ( next_atom[len] == k ) {
                    if ( len < best_len ) {
                        memcpy( best_next_neigh, next_neigh, sizeof(best_next_neigh) );
                        memcpy( best_next_atom,  next_atom,  sizeof(best_next_atom) );
                        best_len = len;
                        num_eql = 0;
                        if ( len == 0 ) {
                            break; /*  path length cannot be smaller than 1 */
                        }
                    } else
                    if ( len == best_len ) {
                        num_eql ++;
                    }
                }
            }
            if ( best_len <= MAX_CUMULENE_LEN && best_next_atom[best_len] == k ) {
                if ( num_eql ) {
                    num_err ++;  /*  program error; no breakpoint here */
                }
                if ( best_len % 2 ) {
                    /*  even number of bonds: chiral atom, draw parity on the cenrtal atom */
                    j = best_next_atom[best_len/2];
                    inf_norm_at[j].cStereoCenterParity = c;
                    str=inf_norm_at[j].at_string;
                    len = strlen(str);
                    if ( len + 3 < (int)sizeof(inf_norm_at[0].at_string) ) {
                        str[len++] = '(';
                        str[len++] = inf_norm_at[j].cStereoCenterParity;
                        str[len++] = ')';
                        str[len] = '\0';
                    }
                } else {
                    /*  odd number of bonds: draw parity on the central bond */
                    if ( best_len == 0 ) {
                        /*  double bond */
                        j = start_at;
                        k = best_next_neigh[0];
                    } else {
                        /*  cumulene */
                        best_len = best_len/2-1;
                        j = best_next_atom[best_len];
                        k = best_next_neigh[best_len+1]; /*  added +1 to display cumulene parity on the middle bond (6-24-2002) */
                    }
                    /*  mark "forward" bond */
                    for ( m = 0; m < MAX_STEREO_BONDS && inf_norm_at[j].cStereoBondParity[m]; m ++ )
                        ;
                    if ( m < MAX_STEREO_BONDS ) {
                        inf_norm_at[j].cStereoBondParity[m] = c;
                        inf_norm_at[j].cStereoBondNumber[m] = k;
                        inf_norm_at[j].cStereoBondWarning[m]   = bAmbiguousStereoBond;
                    } else {
                        num_err ++;  /*  program error; no breakpoint here */
                    }
                    /*  mark "backward" bond */
                    n = norm_at[j].neighbor[k];
                    for ( k = 0; k < norm_at[n].valence && j != (int)norm_at[n].neighbor[k]; k ++ )
                        ;
                    if ( k < norm_at[n].valence ) {
                        j = n;
                        for ( m = 0; m < MAX_STEREO_BONDS && inf_norm_at[j].cStereoBondParity[m]; m ++ )
                            ;
                        if ( m < MAX_STEREO_BONDS ) {
                            inf_norm_at[j].cStereoBondParity[m] = c;
                            inf_norm_at[j].cStereoBondNumber[m] = k;
                            inf_norm_at[j].cStereoBondWarning[m]   = bAmbiguousStereoBond;
                        } else {
                            num_err ++;  /*  program error; no breakpoint here */
                        }
                    } else {
                        num_err ++;  /*  program error; no breakpoint here */
                    }
                }
            } else {
                num_err ++;  /*  program error; no breakpoint here */
            }
        }
    }

    for ( i = 0; i < init_num_at; i ++ ) {
        /*  canonical numbers */
        if ( inf_norm_at[i].nCanonNbr ) {
            str = inf_norm_at[i].at_string;
            len = strlen(str);
            len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nCanonNbr );
            if ( inf_norm_at[i].nCanonEquNbr || inf_norm_at[i].nTautGroupCanonNbr || (inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT) ) {
                if ( inf_norm_at[i].nCanonEquNbr ) {
                    len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nCanonEquNbr );
                } else 
                if ( len + 1 < len_str ) {
                    len += 1;
                    strcat( str, "/" );
                }
            }
            /*  tautomeric groups */
            if ( inf_norm_at[i].nTautGroupCanonNbr ) {
                len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nTautGroupCanonNbr );
                if ( inf_norm_at[i].nTautGroupEquNbr ) {
                    len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nTautGroupEquNbr );
                }
            }
            if ( (inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT) && len+2 <= len_str ) {
                str[len++] = '/';
                str[len++] = '*';
                str[len]   = '\0';
            }
#ifdef DISPLAY_DEBUG_DATA
            if ( inf_norm_at[i].nDebugData ) {
                len += (*MakeNumber)( str+len, len_str-len, "`", (int)inf_norm_at[i].nDebugData );
            }
#endif
        }
    }


exit_function:

    if ( nNormAtNosInCanonOrd )
        inchi_free( nNormAtNosInCanonOrd );


    return ret;
}                              
/***************************************************************************************/
int FillOutOneCanonInfAtom(inp_ATOM *inp_norm_at, INF_ATOM_DATA *inf_norm_at_data,
                           AT_NUMB *pStereoFlags, int init_num_at, int offset, int offset_H, int bIsotopic,
                           INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode)
{
    int          i, j, m, n, num_stereo, k, c, ret, len_str, len, atw, num_err;
    int          next_atom[MAX_CUMULENE_LEN+1], best_next_atom[MAX_CUMULENE_LEN+1], cur_atom;
    int          next_neigh[MAX_CUMULENE_LEN+1], best_next_neigh[MAX_CUMULENE_LEN+1], best_len, bIncludeIsotopicH;
    int          num_iso_H[NUM_H_ISOTOPES];
    char        *str;
    AT_NUMB      g, e;
    int          num_at           = pINChI->nNumberOfAtoms;
    int          num_H            = init_num_at - num_at; /* number of removed H */
    int          nNumberOfTGroups = (pINChI->lenTautomer && pINChI->nTautomer && pINChI->nTautomer[0])? (int)pINChI->nTautomer[0] : 0;
    AT_NUMB     *nOrigAtNosInCanonOrd;
    INChI_Stereo *Stereo;
    AT_NUMB     *nConstitEquNumbers;
    AT_NUMB     *nConstitEquTGroupNumbers;
    S_CHAR      *t_parity = NULL;
    AT_NUMB     *nNumber = NULL;

    AT_NUMB     *nNormAtNosInCanonOrd;
    int (*MakeNumber)(char*, int, const char*, int) = bAbcNumbers? MakeAbcNumber:MakeDecNumber;
    int bRel= (0 != (nMode & ( REQ_MODE_RELATIVE_STEREO)));
    int bRac= (0 != (nMode & ( REQ_MODE_RACEMIC_STEREO)));
    int bRelRac= bRel || bRac;
    int bDoDisplaySp3 = 1;

    inf_ATOM *inf_norm_at   = (inf_norm_at_data && inf_norm_at_data->at)? inf_norm_at_data->at+offset : NULL;
    inp_ATOM *norm_at       = inp_norm_at? inp_norm_at + offset : NULL;
    inf_ATOM *inf_norm_at_H = (inf_norm_at_data && inf_norm_at_data->at)? inf_norm_at_data->at+offset_H : NULL;
    inp_ATOM *norm_at_H     = inp_norm_at? inp_norm_at + offset_H : NULL;

    ret = 0;
    num_err = 0;

    if ( !inf_norm_at )
        return ret;
    /* -- already added in FillOutCompositeCanonInfAtom() --
    if ( bIsotopic ) {
        for ( i = 0, j = 0; i < NUM_H_ISOTOPES; i ++ ) {
            if ( pINChI_Aux->nNumRemovedIsotopicH[i] ) {
                inf_norm_at_data->num_iso_H[i] += pINChI_Aux->nNumRemovedIsotopicH[i];
                inf_norm_at_data->num_removed_iso_H ++;
            }
        }
    }
    */

    if ( bIsotopic && !(pINChI->nNumberOfIsotopicAtoms || pINChI->nNumberOfIsotopicTGroups ||
                        pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0]>1) )
        bIsotopic = 0;
    
    Stereo                   = bIsotopic? pINChI->StereoIsotopic :
                                          pINChI->Stereo;
    bDoDisplaySp3 = (NULL != Stereo) && (Stereo->nNumberOfStereoCenters > 0);

#if( REL_RAC_STEREO_IGN_1_SC == 1 )
    if ( bDoDisplaySp3 && bRelRac && Stereo->nNumberOfStereoCenters < 2 &&
         (Stereo->nCompInv2Abs || ATOM_PARITY_ILL_DEF(Stereo->t_parity[0]) ) ) {
        bDoDisplaySp3 = 0;
        if ( Stereo->nCompInv2Abs ) {
            inf_norm_at_data->StereoFlags |= bRel? INF_STEREO_REL : bRac? INF_STEREO_RAC : 0;
        }
    }
#endif    
    /* flag has stereo */
    if ( (NULL != Stereo) && (bDoDisplaySp3 || Stereo->nNumberOfStereoBonds > 0) ) {
        (*pStereoFlags) |= INF_STEREO;
    }

    /*
    if ( bDoDisplaySp3 && bRelRac && Stereo->nCompInv2Abs && Stereo->nNumberOfStereoCenters < 2 ) {
        bDoDisplaySp3 = 0;
    }
    */
    if ( bDoDisplaySp3 && Stereo->nCompInv2Abs ) {
        /* inversion changes stereo */
        if ( bRel ) {
            (*pStereoFlags) |= INF_STEREO_REL;
        } else
        if ( bRac ) {
            (*pStereoFlags) |= INF_STEREO_RAC;
        } else {
            (*pStereoFlags) |= INF_STEREO_ABS;
        }
        if ( bRelRac ) {
            (*pStereoFlags) |= (Stereo->nCompInv2Abs > 0)? INF_STEREO_NORM : INF_STEREO_INV;
        }
    }
    if ( bDoDisplaySp3 && Stereo->nCompInv2Abs < 0 && !bRelRac ) {
       /* display Inv stereo which is Absolute Stereo */
        nNumber  = Stereo->nNumberInv;
        t_parity = Stereo->t_parityInv;
        nOrigAtNosInCanonOrd     = bIsotopic? pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv :
                                              pINChI_Aux->nOrigAtNosInCanonOrdInv;
    } else {
       /* display Output stereo which is Absolute Stereo */
        if ( bDoDisplaySp3 ) {
            nNumber  = Stereo->nNumber;
            t_parity = Stereo->t_parity;
        }
        nOrigAtNosInCanonOrd     = bIsotopic? pINChI_Aux->nIsotopicOrigAtNosInCanonOrd :
                                              pINChI_Aux->nOrigAtNosInCanonOrd;
    }

    nConstitEquNumbers       = bIsotopic? pINChI_Aux->nConstitEquIsotopicNumbers :
                                          pINChI_Aux->nConstitEquNumbers;
    nConstitEquTGroupNumbers = bIsotopic? pINChI_Aux->nConstitEquIsotopicTGroupNumbers :
                                          pINChI_Aux->nConstitEquTGroupNumbers;
    memset( inf_norm_at,   0, num_at*sizeof(inf_norm_at[0]) );
    if ( num_H > 0 ) {
        memset( inf_norm_at_H, 0, num_H*sizeof(inf_norm_at[0]) );
    }

    /*  obtain norm_at[] atom numbers (from zero) in order of canonical numbers */
    nNormAtNosInCanonOrd = (AT_NUMB *)inchi_calloc( num_at, sizeof(nNormAtNosInCanonOrd[0]) );
    if ( ret = GetAtomOrdNbrInCanonOrd( norm_at, nNormAtNosInCanonOrd, nOrigAtNosInCanonOrd, num_at ) ) {
        goto exit_function;
    }
    
    /*  atom canonical and equivalence numbers > 0 */
    for ( i = 0; i < num_at; i ++ ) {
        j = (int)nNormAtNosInCanonOrd[i];
        if ( j < 0 || j >= num_at )
            continue;
        inf_norm_at[j].nCanonNbr = (AT_NUMB)(i+1);
        inf_norm_at[j].nCanonEquNbr = nConstitEquNumbers[i];
#ifdef DISPLAY_DEBUG_DATA
        inf_norm_at[j].nDebugData = 0;
#if( DISPLAY_DEBUG_DATA == DISPLAY_DEBUG_DATA_C_POINT )
        inf_norm_at[j].nDebugData = norm_at[j].c_point;
#endif
#endif
    }
    /*  tautomeric groups */
    if ( nNumberOfTGroups ) {
        /*
         :   start from 1: bypass number of t-groups
         :   j is a counter within the current t-group
         :   g is a tautomeric group canonical number
         :   e is a tautomeric group equivalence
         */
        for ( g = 1, i = 1; g <= nNumberOfTGroups; g ++ ) {
            n = (int)pINChI->nTautomer[i] - INCHI_T_NUM_MOVABLE; /*  number of atoms in t-group */
            e = nConstitEquTGroupNumbers[(int)g - 1];
            /*  bypass number of hydrogen atoms, negative charges, ... */
            for ( i += INCHI_T_NUM_MOVABLE+1, j = 0; j < n && i < pINChI->lenTautomer; j ++, i ++ ) {
                /*  scan canonical numbers of atoms within the atom t-group */
                k = (int)nNormAtNosInCanonOrd[(int)pINChI->nTautomer[i]-1];
                inf_norm_at[k].nTautGroupCanonNbr = g;
                inf_norm_at[k].nTautGroupEquNbr   = e;
            }
        }
        if ( i != pINChI->lenTautomer || g != nNumberOfTGroups+1 ) {
            ret = CT_TAUCOUNT_ERR;  /*   <BRKPT> */
            goto exit_function;
        }
    }
    /* atoms that may exchange isotopic H */
    if ( bIsotopic && pINChI->nPossibleLocationsOfIsotopicH && (n = (int)pINChI->nPossibleLocationsOfIsotopicH[0]) ) {
        for ( i = 1; i < n; i ++ ) {
            j = (int)pINChI->nPossibleLocationsOfIsotopicH[i];
            k = (int)nNormAtNosInCanonOrd[j - 1];
            if ( !inf_norm_at[k].nTautGroupCanonNbr ) {
                inf_norm_at[k].cFlags |= AT_FLAG_ISO_H_POINT;
            }
        }
    }
#if( DISPLAY_RING_SYSTEMS == 1 )
    /*  debug only */
    for ( j = 0; j < num_at; j ++ ) {
        inf_norm_at[j].nCanonNbr          = norm_at[j].nBlockSystem;
        inf_norm_at[j].nCanonEquNbr       = norm_at[j].nRingSystem;
#if( USE_DISTANCES_FOR_RANKING == 1 )
        inf_norm_at[j].nTautGroupCanonNbr   = norm_at[j].nDistanceFromTerminal;
        inf_norm_at[j].nTautGroupEquNbr     = norm_at[j].bCutVertex;
#else
        inf_norm_at[j].nTautGroupCanonNbr   = norm_at[j].bCutVertex;
        inf_norm_at[j].nTautGroupEquNbr   = 0;
#endif
    }
#endif



    /*  Write isotopic mass, chemical element symbols and hydrogens, charge, radical, canon. numbers */
    len_str = sizeof(inf_norm_at[0].at_string);
    for ( i = 0; i < init_num_at; i ++ ) {
        inf_ATOM *cur_inf_norm_at = (i < num_at)? inf_norm_at+i : inf_norm_at_H+i-num_at;
        inp_ATOM *cur_norm_at     = (i < num_at)? norm_at    +i : norm_at_H    +i-num_at;
        str = cur_inf_norm_at->at_string;
        len = 0;
        bIncludeIsotopicH = bIsotopic && (i >= num_at || !inf_norm_at[i].nTautGroupCanonNbr && !(inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT));
        /*  isotopic mass */
        atw = 0;
        if ( cur_norm_at->iso_atw_diff && bIsotopic ) {
            if ( cur_norm_at->at_type == ATT_PROTON ) {
                ; /* do not set isotopic mass of a tautomeric proton */
            } else
            if ( num_at <= i && cur_norm_at->el_number == PERIODIC_NUMBER_H && cur_norm_at->chem_bonds_valence == 1 &&
                 !cur_norm_at->charge && !cur_norm_at->radical && !cur_norm_at->num_H &&
                 0 <= (j=(int)cur_norm_at->neighbor[0]-offset) && j < num_at &&
                 (inf_norm_at[j].nTautGroupCanonNbr || (inf_norm_at[j].cFlags & AT_FLAG_ISO_H_POINT) ) ) {
                ; /* do not set isotopic mass of an exchangeable proton */
            } else {
                atw = get_atw(cur_norm_at->elname);
                atw += (cur_norm_at->iso_atw_diff>0)? cur_norm_at->iso_atw_diff-1:cur_norm_at->iso_atw_diff;
                /*len += sprintf( str+len, "^%d", atw );*/
            }
        }
        /*  element name */
        if ( cur_norm_at->el_number == PERIODIC_NUMBER_H && 2 <= atw && atw <= 3 ) {
            len += sprintf( str+len, "%s", atw==2? "D" : "T" );
        } else {
            if ( atw ) {
                len += sprintf( str+len, "^%d", atw );
            }
            len += sprintf( str+len, "%s", cur_norm_at->elname );
        }
        /*  hydrogens */
        /*  find number of previuosly removed terminal hydrogen atoms */
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
            num_iso_H[j] = cur_norm_at->num_iso_H[j];
        }
        for ( j = 0, n = (int)cur_norm_at->num_H; j < num_H; j ++ ) {
            /*  subtract number of removed terminal */
            /*  H atoms from the total number of H atoms */
            if ( i == (int)norm_at_H[j].neighbor[0]-offset ) {
                n -= 1;
                m = (int)norm_at_H[j].iso_atw_diff-1;
                if ( 0 <= m && m < NUM_H_ISOTOPES ) {
                    /*  subtract number of removed terminal isotopic */
                    /*  H atoms from the total number of isotopic H atoms */
                    num_iso_H[m] -= 1;
                }
            }
        }
        if ( bIncludeIsotopicH ) {
            /*  subtract number of isotopic H atoms from the total number of H atoms */
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                n -= num_iso_H[j];
            }
        }
        /*  non-isotopic hydrogen atoms */
        if ( n > 1 ) {
            len += sprintf( str+len, "H%d", n );
        } else
        if ( n == 1 ) {
            len += sprintf( str+len, "H" );
        }
        /*  isotopic hydrogen atoms */
        if ( bIncludeIsotopicH ) {
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                if ( num_iso_H[j] ) {
                    if ( j == 0 || j != 1 && j != 2 ) {
                        len += sprintf( str+len, "^%dH", j+1 );
                    } else {
                        len += sprintf( str+len, j == 1? "D" : "T" );
                    }
                    if ( num_iso_H[j] != 1 ) {
                        len += sprintf( str+len, "%d", (int)num_iso_H[j] );
                    }
                }
            }
        }
        if ( cur_norm_at->el_number == PERIODIC_NUMBER_H && str[0] == str[1] ) {
            char *q;
            if ( !str[2] ) {
                str[1] = '2';  /* quick fix: replace HH with H2 */
            } else
            if ( isdigit( UCINT str[2] ) && (n = strtol( str+2, &q, 10 )) && !q[0] ) {
                len = 1 + sprintf( str+1, "%d", n+1 );
            }
        }
        /*  charge */
        if ( abs(cur_norm_at->charge) > 1 )
            len += sprintf( str+len, "%+d", cur_norm_at->charge );
        else
        if ( abs(cur_norm_at->charge) == 1 )
            len += sprintf( str+len, "%s", cur_norm_at->charge>0? "+" : "-" );
        /*  radical */
        if ( cur_norm_at->radical )
            len += sprintf( str+len, "%s", cur_norm_at->radical==RADICAL_SINGLET? ":" :
                                           cur_norm_at->radical==RADICAL_DOUBLET? "." :
                                           cur_norm_at->radical==RADICAL_TRIPLET? ".." : "?");
    }

    /*  stereogenic centers */
    if ( bDoDisplaySp3 && Stereo && 0 < (num_stereo = Stereo->nNumberOfStereoCenters) ) {
        for ( i = 0; i < num_stereo; i ++ ) {
            j = (int)nNormAtNosInCanonOrd[(int)nNumber[i]-1];
            c = t_parity[i];
            c = c==1? '-' : c==2? '+' : c==3? 'u' : c== 4? '?' : '*';
            inf_norm_at[j].cStereoCenterParity = c;
            str=inf_norm_at[j].at_string;
            len = strlen(str);
            if ( len + 3 < (int)sizeof(inf_norm_at[0].at_string) ) {
                str[len++] = '(';
                str[len++] = inf_norm_at[j].cStereoCenterParity;
                str[len++] = ')';
                str[len] = '\0';
                /*  mark ambuguous stereo center */
                if ( norm_at[j].bAmbiguousStereo && (c=='+' || c=='-' || c=='?') &&  str[0] != '!' &&
                     len+1 < (int)sizeof(inf_norm_at[0].at_string) ) {
                    memmove( str+1, str, len+1 );
                    str[0] = '!'; /* output the atom in red color */
                }
            }
        }
    }
    
    /*  stereogenic bonds */
    /*  (cumulenes with odd number of double bonds are stereocenters, */
    /*   and atom parity should be set) */
    if ( Stereo && 0 < (num_stereo = Stereo->nNumberOfStereoBonds) ) {
        for ( i = 0; i < num_stereo; i ++ ) {
            int start_at, num_eql=0, bAmbiguousStereoBond = 0;
            j = (int)nNormAtNosInCanonOrd[(int)Stereo->nBondAtom1[i]-1];
            k = (int)nNormAtNosInCanonOrd[(int)Stereo->nBondAtom2[i]-1];
            start_at = j;
            c = Stereo->b_parity[i];
                                
            c = c==1? '-' : c==2? '+' : c==3? 'u' : c== 4? '?' : '*';

            /*  mark ambuguous stereo bond atom(s) */
            if ( norm_at[j].bAmbiguousStereo && (c=='+' || c=='-' ) &&
                 (len=strlen(str=inf_norm_at[j].at_string)+1) < (int)sizeof(inf_norm_at[0].at_string) &&
                  str[0] != '!' ) {
                memmove( str+1, str, len );
                str[0] = '!'; /* output the atom in red color */
                bAmbiguousStereoBond ++;
            }
            if ( norm_at[k].bAmbiguousStereo && (c=='+' || c=='-') &&
                 (len=strlen(str=inf_norm_at[k].at_string)+1) < (int)sizeof(inf_norm_at[0].at_string) &&
                  str[0] != '!' ) {
                memmove( str+1, str, len );
                str[0] = '!'; /* output the atom in red color */
                bAmbiguousStereoBond ++;
            }

            /*  find the opposite atom k. */
            /*  Note: since it may be a cumulene, find the shortest(best) path */
            /*  to atom number k to avoid confusion in case of, for example, */
            /*  4-member aromatic rings. */
            best_len = MAX_CUMULENE_LEN+1; /* moved here from inside the cycle 1-8-2003 */
            for ( n = 0; n < norm_at[j].valence; n ++ ) {
                if ( norm_at[j].bond_type[n] == BOND_SINGLE ) {
                    /*  single bond cannot be stereogenic. */
                    continue;
                }
                /* best_len = MAX_CUMULENE_LEN+1; */
                len = 0; /*  number of bonds in cumulene - 1 */
                next_atom[len]  = (int)norm_at[j].neighbor[n]-offset;
                next_neigh[len] = n;
                cur_atom = j;
                while ( next_atom[len] != k && len < MAX_CUMULENE_LEN && 2 == norm_at[next_atom[len]].valence ) {
                    next_neigh[len+1] = ((int)norm_at[next_atom[len]].neighbor[0]-offset == cur_atom);
                    next_atom[len+1]  =  (int)norm_at[next_atom[len]].neighbor[next_neigh[len+1]]-offset;
                    cur_atom = next_atom[len];
                    len ++;
                }
                if ( next_atom[len] == k ) {
                    if ( len < best_len ) {
                        memcpy( best_next_neigh, next_neigh, sizeof(best_next_neigh) );
                        memcpy( best_next_atom,  next_atom,  sizeof(best_next_atom) );
                        best_len = len;
                        num_eql = 0;
                        if ( len == 0 ) {
                            break; /*  path length cannot be smaller than 1 */
                        }
                    } else
                    if ( len == best_len ) {
                        num_eql ++;
                    }
                }
            }
            if ( best_len <= MAX_CUMULENE_LEN && best_next_atom[best_len] == k ) {
                if ( num_eql ) {
                    num_err ++;  /*  program error; no breakpoint here */
                }
                if ( best_len % 2 ) {
                    /*  even number of bonds: chiral atom, draw parity on the cenrtal atom */
                    j = best_next_atom[best_len/2];
                    inf_norm_at[j].cStereoCenterParity = c;
                    str=inf_norm_at[j].at_string;
                    len = strlen(str);
                    if ( len + 3 < (int)sizeof(inf_norm_at[0].at_string) ) {
                        str[len++] = '(';
                        str[len++] = inf_norm_at[j].cStereoCenterParity;
                        str[len++] = ')';
                        str[len] = '\0';
                    }
                } else {
                    /*  odd number of bonds: draw parity on the central bond */
                    if ( best_len == 0 ) {
                        /*  double bond */
                        j = start_at;
                        k = best_next_neigh[0];
                    } else {
                        /*  cumulene */
                        best_len = best_len/2-1;
                        j = best_next_atom[best_len];
                        k = best_next_neigh[best_len+1]; /*  added +1 to display cumulene parity on the middle bond (6-24-2002) */
                    }
                    /*  mark "forward" bond */
                    for ( m = 0; m < MAX_STEREO_BONDS && inf_norm_at[j].cStereoBondParity[m]; m ++ )
                        ;
                    if ( m < MAX_STEREO_BONDS ) {
                        inf_norm_at[j].cStereoBondParity[m] = c;
                        inf_norm_at[j].cStereoBondNumber[m] = k;
                        inf_norm_at[j].cStereoBondWarning[m]  = bAmbiguousStereoBond;
                    } else {
                        num_err ++;  /*  program error; no breakpoint here */
                    }
                    /*  mark "backward" bond */
                    n = norm_at[j].neighbor[k]-offset;
                    for ( k = 0; k < norm_at[n].valence && j != (int)norm_at[n].neighbor[k]-offset; k ++ )
                        ;
                    if ( k < norm_at[n].valence ) {
                        j = n;
                        for ( m = 0; m < MAX_STEREO_BONDS && inf_norm_at[j].cStereoBondParity[m]; m ++ )
                            ;
                        if ( m < MAX_STEREO_BONDS ) {
                            inf_norm_at[j].cStereoBondParity[m] = c;
                            inf_norm_at[j].cStereoBondNumber[m] = k;
                            inf_norm_at[j].cStereoBondWarning[m]   = bAmbiguousStereoBond;
                        } else {
                            num_err ++;  /*  program error; no breakpoint here */
                        }
                    } else {
                        num_err ++;  /*  program error; no breakpoint here */
                    }
                }
            } else {
                num_err ++;  /*  program error; no breakpoint here */
            }
        }
    }

    for ( i = 0; i < num_at; i ++ ) {
        /*  canonical numbers */
        if ( inf_norm_at[i].nCanonNbr ) {
            str = inf_norm_at[i].at_string;
            len = strlen(str);
            len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nCanonNbr );
            if ( inf_norm_at[i].nCanonEquNbr || inf_norm_at[i].nTautGroupCanonNbr || (inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT) ) {
                if ( inf_norm_at[i].nCanonEquNbr ) {
                    len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nCanonEquNbr );
                } else 
                if ( len + 1 < len_str ) {
                    len += 1;
                    strcat( str, "/" );
                }
            }
            /*  tautomeric groups */
            if ( inf_norm_at[i].nTautGroupCanonNbr ) {
                len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nTautGroupCanonNbr );
                if ( inf_norm_at[i].nTautGroupEquNbr ) {
                    len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_norm_at[i].nTautGroupEquNbr );
                }
            }
            if ( (inf_norm_at[i].cFlags & AT_FLAG_ISO_H_POINT) && len+2 <= len_str ) {
                str[len++] = '/';
                str[len++] = '*';
                str[len]   = '\0';
            }
#ifdef DISPLAY_DEBUG_DATA
            if ( inf_norm_at[i].nDebugData ) {
                len += (*MakeNumber)( str+len, len_str-len, "`", (int)inf_norm_at[i].nDebugData );
            }
#endif
        }
    }


exit_function:

    if ( nNormAtNosInCanonOrd )
        inchi_free( nNormAtNosInCanonOrd );

    return ret;
}

/***************************************************************************************/
int FillOutInputInfAtom(inp_ATOM *inp_at, INF_ATOM_DATA *inf_at_data, int init_num_at, int num_removed_H,
                        int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, int bAbcNumbers)
{
    int          i, j, m, n, ret, len_str, len, atw;
    int          num_iso_H[NUM_H_ISOTOPES];
    char        *str;
    int          num_at           = init_num_at - num_removed_H;
    int (*MakeNumber)(char*, int, const char*, int) = MakeDecNumber;

    inf_ATOM *inf_at = inf_at_data? inf_at_data->at : NULL;


    ret = 0;
    

    if ( !inf_at )
        return ret;

    memset( inf_at, 0, init_num_at*sizeof(inf_at[0]) );

    inf_at_data->nNumRemovedProtons = nNumRemovedProtons;
    MakeRemovedProtonsString( nNumRemovedProtons, nNumRemovedProtonsIsotopic, NULL, bIsotopic, inf_at_data->szRemovedProtons, NULL );
    /*  atom canonical and equivalence numbers > 0 */
    for ( i = 0; i < num_at; i ++ ) {
#if( DISPLAY_ORIG_AT_NUMBERS == 1 )
        inf_at[i].nCanonNbr = inp_at[i].orig_at_number;
#else
        inf_at[i].nCanonNbr = (AT_NUMB)(i+1);
#endif
    }
    /*  Write isotopic mass, chemical element symbols and hydrogens, charge, radical, canon. numbers */
    len_str = sizeof(inf_at[0].at_string);
    for ( i = 0; i < init_num_at; i ++ ) {
        str = inf_at[i].at_string;
        len = 0;
        /*  isotopic mass */
        atw = 0;
        if ( inp_at[i].iso_atw_diff && bIsotopic ) {
            atw = get_atw(inp_at[i].elname);
            atw += (inp_at[i].iso_atw_diff>0)? inp_at[i].iso_atw_diff-1:inp_at[i].iso_atw_diff;
            /*len += sprintf( str+len, "^%d", atw );*/
        }
        /*  element name */
        if ( inp_at[i].el_number == PERIODIC_NUMBER_H && 2 <= atw && atw <= 3 ) {
            len += sprintf( str+len, "%s", atw==2? "D" : "T" );
        } else {
            if ( atw ) {
                len += sprintf( str+len, "^%d", atw );
            }
            len += sprintf( str+len, "%s", inp_at[i].elname );
        }
        /*  hydrogens */
        /*  find number of previuosly removed terminal hydrogen atoms */
        for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
            num_iso_H[j] = inp_at[i].num_iso_H[j];
        }
        for ( j = num_at, n = (int)inp_at[i].num_H; j < init_num_at; j ++ ) {
            /*  subtract number of removed terminal */
            /*  H atoms from the total number of H atoms */
            if ( i == (int)inp_at[j].neighbor[0] ) {
                n -= 1;
                m = (int)inp_at[j].iso_atw_diff-1;
                if ( 0 <= m && m < NUM_H_ISOTOPES ) {
                    /*  subtract number of removed terminal isotopic */
                    /*  H atoms from the total number of isotopic H atoms */
                    num_iso_H[m] -= 1;
                }
            }
        }
        if ( bIsotopic && !bAdd_DT_to_num_H ) {
            /*  subtract number of isotopic H atoms from the total number of H atoms */
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                n -= num_iso_H[j];
            }
        }
        /*  non-isotopic hydrogen atoms */
        if ( n > 1 ) {
            len += sprintf( str+len, "H%d", n );
        } else
        if ( n == 1 ) {
            len += sprintf( str+len, "H" ); /* fixed 12-21-2002: removed 3rd argument */
        }
        if ( bIsotopic ) {
            /*  isotopic hydrogen atoms */
            for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                if ( num_iso_H[j] ) {
                    if ( j == 0 || j != 1 && j != 2 ) {
                        len += sprintf( str+len, "^%dH", j+1 );
                    } else {
                        len += sprintf( str+len, j == 1? "D" : "T" );
                    }
                    if ( num_iso_H[j] != 1 ) {
                        len += sprintf( str+len, "%d", (int)num_iso_H[j] );
                    }
                }
            }
        }
        if ( inp_at[i].el_number == PERIODIC_NUMBER_H && str[0] == str[1] ) {
            char *q;
            if ( !str[2] ) {
                str[1] = '2';  /* quick fix: replace HH with H2 */
            } else
            if ( isdigit( UCINT str[2] ) && (n = strtol( str+2, &q, 10 )) && !q[0] ) {
                len = 1 + sprintf( str+1, "%d", n+1 );
            }
        }
        /*
        if (  str[0] == 'H' && str[1] == 'H' && !str[2] ) {
            str[1] = '2';
        }
        */
        /*  charge */
        if ( abs(inp_at[i].charge) > 1 )
            len += sprintf( str+len, "%+d", inp_at[i].charge );
        else
        if ( abs(inp_at[i].charge) == 1 )
            len += sprintf( str+len, "%s", inp_at[i].charge>0? "+" : "-" );
        /*  radical */
        if ( inp_at[i].radical )
            len += sprintf( str+len, "%s", inp_at[i].radical==RADICAL_SINGLET? ":" :
                                           inp_at[i].radical==RADICAL_DOUBLET? "." :
                                           inp_at[i].radical==RADICAL_TRIPLET? ".." : "?");
    }

    for ( i = 0; i < init_num_at; i ++ ) {
        /*  canonical numbers */
        if ( inf_at[i].nCanonNbr ) {
            str = inf_at[i].at_string;
            len = strlen(str);
            len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_at[i].nCanonNbr );
            if ( inf_at[i].nCanonEquNbr || inf_at[i].nTautGroupCanonNbr ) {
                len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_at[i].nCanonEquNbr );
            }
            /*  tautomeric groups */
            if ( inf_at[i].nTautGroupCanonNbr ) {
                len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_at[i].nTautGroupCanonNbr );
                if ( inf_at[i].nTautGroupEquNbr ) {
                    len += (*MakeNumber)( str+len, len_str-len, "/", (int)inf_at[i].nTautGroupEquNbr );
                }
            }
        }
    }
    ret = init_num_at;

    return ret;
}
/**********************************************************************************************/
int FillOutInfAtom(inp_ATOM *norm_at, INF_ATOM_DATA *inf_norm_at_data, int init_num_at, int num_removed_H,
                   int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic,
                   INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode )
{
    if ( norm_at && inf_norm_at_data && inf_norm_at_data->at ) {
        if ( pINChI && pINChI_Aux ) {
            return FillOutCanonInfAtom( norm_at, inf_norm_at_data, init_num_at, bIsotopic, pINChI,
                                        pINChI_Aux, bAbcNumbers, nMode);
        } else {
            return FillOutInputInfAtom( norm_at, inf_norm_at_data, init_num_at, num_removed_H, bAdd_DT_to_num_H,
                                        nNumRemovedProtons, nNumRemovedProtonsIsotopic, bIsotopic, bAbcNumbers);
        }
    }
    return 0;
}
/***************************************************************************************/
int FillOutCompositeCanonInfAtom(COMP_ATOM_DATA *composite_norm_data, INF_ATOM_DATA *inf_norm_at_data,
                                 int bIsotopic, int bTautomeric,
                                 PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2, int bAbcNumbers, INCHI_MODE nMode)
{
    int i, num_components, j, k, ret;
    inp_ATOM *inp_norm_at;
    INChI     *pINChI;
    INChI_Aux *pINChI_Aux;
    int      num_inp_at, num_at, num_H, offset, offset_H, next_offset, next_offset_H;

    if ( composite_norm_data && inf_norm_at_data && (bTautomeric == TAUT_INI || pINChI2 && pINChI_Aux2) ) {
        composite_norm_data += bTautomeric;
        inp_norm_at          = composite_norm_data->at;
        num_components       = composite_norm_data->num_components;
        offset   = 0;
        offset_H = composite_norm_data->num_at - composite_norm_data->num_removed_H;
        if ( bTautomeric == TAUT_INI ) {
            ret = FillOutInputInfAtom( composite_norm_data->at, inf_norm_at_data, composite_norm_data->num_at,
                                       composite_norm_data->num_removed_H, 0 /*bAdd_DT_to_num_H*/,
                                       composite_norm_data->nNumRemovedProtons,
                                       composite_norm_data->nNumRemovedProtonsIsotopic,
                                       bIsotopic, bAbcNumbers);
            return ret;
        } else {
            for ( i = 0; i < num_components; i ++ ) {
                j = inchi_min(bTautomeric, TAUT_YES);
                /* count isotopic H on removed atoms -- isolated H(+) cations */
                inf_norm_at_data->nNumRemovedProtons += pINChI_Aux2[i][j]->nNumRemovedProtons;
                if ( bIsotopic && bTautomeric == TAUT_YES ) {
                    for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                        if ( pINChI_Aux2[i][j]->nNumRemovedIsotopicH[k] ) {
                            inf_norm_at_data->num_iso_H[k] += pINChI_Aux2[i][j]->nNumRemovedIsotopicH[k];
                            inf_norm_at_data->num_removed_iso_H += pINChI_Aux2[i][j]->nNumRemovedIsotopicH[k];
                        }
                    }
                }
                /* ignore deleted components */
                if ( pINChI2[i][j] && pINChI2[i][j]->bDeleted ) {
                    continue;
                }
                if ( !pINChI2[i][j] || !pINChI2[i][j]->nNumberOfAtoms ) {
                    j = ALT_TAUT(j);
                    if ( !pINChI2[i][j] || !pINChI2[i][j]->nNumberOfAtoms ) {
                        continue; /* error ??? */
                    }
                }
                pINChI     = pINChI2[i][j];
                pINChI_Aux = pINChI_Aux2[i][j];
                next_offset   = composite_norm_data->nOffsetAtAndH[2*i];
                next_offset_H = composite_norm_data->nOffsetAtAndH[2*i+1];
                num_at   = next_offset - offset;
                if ( num_at <= 0 )
                    continue;
                num_H    = next_offset_H - offset_H;
                num_inp_at = num_at + num_H;
                if ( num_at != pINChI->nNumberOfAtoms || num_at != pINChI_Aux->nNumberOfAtoms ) {
                    return 0; /* error */
                }
                ret =  FillOutOneCanonInfAtom(inp_norm_at, inf_norm_at_data,
                               inf_norm_at_data->pStereoFlags+i+1, num_inp_at,
                               offset, offset_H, bIsotopic, pINChI, pINChI_Aux, bAbcNumbers, nMode);
                if ( ret )
                    return 0; /* error */

                inf_norm_at_data->StereoFlags |= inf_norm_at_data->pStereoFlags[i+1];
                offset   = next_offset;
                offset_H = next_offset_H;
            }
        }
        MakeRemovedProtonsString( inf_norm_at_data->nNumRemovedProtons, inf_norm_at_data->num_iso_H, NULL, bIsotopic,
                                  inf_norm_at_data->szRemovedProtons, &inf_norm_at_data->num_removed_iso_H );
    }
    return 1;
}
#endif /* } ifndef INCHI_ANSI_ONLY */
/**********************************************************************************************/
int CheckCanonNumberingCorrectness(int num_atoms, int num_at_tg,
                 sp_ATOM *at, CANON_STAT *pCS, int bTautomeric,
                 char *pStrErrStruct )
{
    int i, ret=0;
    AT_NUMB *pCanonOrd=NULL;
    int nErrorCode = 0;
    AT_NUMB *pCanonRank; /* canonical ranks of the atoms or tautomeric groups */
    AT_NUMB *pCanonRankAtoms=NULL;
    
    static int count=0; /* for debug only */
    count ++;

    pCanonRankAtoms = (AT_NUMB *)inchi_calloc( num_at_tg+1, sizeof(pCanonRankAtoms[0]) );
    
    /**********************************************************************************************
     *
     *  non-isotopic part
     */
    pCanonOrd            = pCS->nLenCanonOrdStereo > 0? pCS->nCanonOrdStereo :
                           pCS->nLenCanonOrd       > 0? pCS->nCanonOrd : NULL;
    pCanonRank           = pCanonRankAtoms;
    if ( pCanonOrd && pCanonRank ) {
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRank[pCanonOrd[i]] = (AT_NUMB)(i+1);
        }
        ret = UpdateFullLinearCT( num_atoms, num_at_tg, at, pCanonRank, pCanonOrd, pCS, 0 );
        if ( ret /*|| memcmp(pCS->LinearCT, pCS->LinearCT2, sizeof(AT_RANK) * pCS->nLenLinearCT )*/ ) {
            nErrorCode     |= WARN_FAILED_STEREO;
        }

    } else {
        nErrorCode  |= ERR_NO_CANON_RESULTS;
        goto exit_function;
    }
    /**********************************************************************************************
     *
     *  isotopic part
     */
    pCanonOrd   = pCS->nLenCanonOrdIsotopicStereo > 0? pCS->nCanonOrdIsotopicStereo :
                  pCS->nLenCanonOrdIsotopic       > 0? pCS->nCanonOrdIsotopic : NULL;
    pCanonRank  = pCanonRankAtoms;

    if ( pCanonOrd && pCanonRank ) {
        for ( i = 0; i < num_at_tg; i ++ ) {
            pCanonRank[pCanonOrd[i]] = (AT_NUMB)(i+1);
        }
        ret = UpdateFullLinearCT( num_atoms, num_at_tg, at, pCanonRank, pCanonOrd, pCS, 0 );
        if ( ret /*|| memcmp(pCS->LinearCT, pCS->LinearCT2, sizeof(AT_RANK) * pCS->nLenLinearCT )*/ ) {
            nErrorCode     |= (pCS->nLenCanonOrdIsotopicStereo? WARN_FAILED_ISOTOPIC_STEREO : WARN_FAILED_ISOTOPIC);
        }

    }

exit_function:
    if ( pCanonRankAtoms )
        inchi_free( pCanonRankAtoms );

    if ( nErrorCode ) {
        return CT_CANON_ERR;
    }
    return 0;
}
